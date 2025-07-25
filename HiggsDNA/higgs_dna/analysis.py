import os
import time
import datetime
import copy
import json
import pickle
import dill
import numpy

import uproot
import awkward

import logging
import orjson
# logger = logging.getLogger(__name__)
from higgs_dna.utils.logger_utils import simple_logger
logger = simple_logger(__name__)

from higgs_dna.samples.sample_manager import SampleManager
from higgs_dna.samples.sample import Sample
from higgs_dna.job_management.managers import LocalManager, CondorManager
from higgs_dna.job_management.task import Task
from higgs_dna.systematics.systematics_producer import SystematicsProducer
from higgs_dna.taggers.tag_sequence import TagSequence
from higgs_dna.utils.misc_utils import load_config, update_dict, is_json_serializable
from higgs_dna.constants import NOMINAL_TAG, CENTRAL_WEIGHT, BRANCHES
from higgs_dna.utils.metis_utils import do_cmd
from higgs_dna.taggers.duplicated_samples_tagger import DuplicatedSamplesTagger
from higgs_dna.taggers.mc_overlap_tagger import MCOverlapTagger
from higgs_dna.systematics.photon_systematics import photon_scale_smear_run3
from higgs_dna.systematics.lepton_systematics import electron_scale_smear_run3

condor=False
def run_analysis(config):
    """
    Function that gets run for each individual job. Performs the following:
        1. Load events from input files
        2. Add any relevant ``Sample`` metadata to events
        3. Produce all ``WeightSystematic``s and ``SystematicWithIndependentCollection``s that **do not** rely on selections performed by or fields added by ``Tagger``s.
        4. Apply ``TagSequence``
        5. Produce all remaining ``WeightSystematic``s and ``SystematicWithIndependentCollection``s that **do** rely on selections performed by or fields added by ``Tagger``s.
        6. Write output events 

    In general, you should not have to call this function yourself.
    It will be configured by the ``AnalysisManager`` for each job.

    :param config: dictionary/json that specifies the physics content (samples, tag sequence, systematics) of a given job
    :type config: dict or str
    """
    t_start = time.time()

    config = load_config(config)
    job_summary = { 
            "config" : config
    }

    ### 1. Load events ###
    t_start_load = time.time()
    events, sum_weights = AnalysisManager.load_events(config)

    # Optional branch mapping in case you have different naming schemes, e.g. you want MET_T1smear_pt to be recast as MET_pt
    # Can be separate for data and MC
    if "branch_map" in config.keys():
        if config["sample"]["is_data"]:
            branch_map = config["branch_map"]["data"]
        else:
            branch_map = config["branch_map"]["mc"]

        if branch_map:
            for x in branch_map:
                logger.debug("[run_analysis] Replacing %s with %s." % (str(x[0]), str(x[1])))
                print("[DEBUG]: ", str(x[0]), list, tuple(x[0]), tuple(x[1]))
                print("[DEBUG]: ", awkward.fields(events))
                print("[DEBUG]: ", awkward.fields(events["Jet"]))
                if isinstance(x[0], list):
                    events[tuple(x[0])] = events[tuple(x[1])]
                else:
                    events[x[0]] = events[x[1]]

    else:
        logger.debug("[run_analysis] No branch map.")
    if not os.path.exists(config["dir"]):
        output_dir = ""
        config["summary_file"] = os.path.split(config["summary_file"])[-1]
    else:
        output_dir = os.path.abspath(config["dir"]) + "/"

    # Record n_events and sum_weights for scale1fb calculation
    output_name = output_dir + config["output_name"]
    job_summary["n_events"] = len(events)
    job_summary["sum_weights"] = sum_weights
    t_elapsed_load = time.time() - t_start_load

    ### 2. Add relevant sample metadata to events ###
    t_start_samples = time.time()
    sample = Sample(**config["sample"])
    events = sample.prep(events)
    t_elapsed_samples = time.time() - t_start_samples

    ### 3. Produce systematics ###
    t_start_syst = time.time()
    systematics_producer = SystematicsProducer(
        name = config["name"],
        options = config["systematics"],
        sample = sample
    )
    events = systematics_producer.produce_weights(events)
    tag_sequence = TagSequence(
        name = config["name"],
        tag_list = config["tag_sequence"],
        sample = sample,
        output_dir=config["output_dir"])
    job_summary["outputs"] = {}
    job_summary["n_events_selected"] = {}
    t_elapsed_syst = 0.
    t_elapsed_taggers = 0.

    # Systematics variations
    for syst_name, ic_syst in systematics_producer.independent_collections.items():
        if not systematics_producer.do_variations:
            continue
        current_time = time.time()
        ics = ic_syst.produce(events)
        t_elapsed_syst += time.time() - current_time        

        for ic_name, events_ic in ics.items():
            # If the IC modifies the nominal value of a branch, update this in the nominal events
            if ic_name == "nominal":
                events = events_ic
                continue

            # Otherwise, this is an up/down variation
            syst_tag = syst_name + "_" + ic_name

            current_time = time.time()

            events_ic, tag_idx_map = tag_sequence.run(events_ic, syst_tag)
            t_elapsed_taggers += time.time() - current_time

            current_time = time.time()
            events_ic = systematics_producer.apply_remaining_weight_systs(events_ic, syst_tag, tag_idx_map)
            t_elapsed_syst += time.time() - current_time

            job_summary["outputs"][syst_tag] = AnalysisManager.write_events(events_ic, config["variables_of_interest"], output_name, syst_tag)       
            job_summary["n_events_selected"][syst_tag] = len(events_ic) 

    # Nominal events
    events, tag_idx_map = tag_sequence.run(events, NOMINAL_TAG)
    events = systematics_producer.apply_remaining_weight_systs(events, NOMINAL_TAG, tag_idx_map)
    job_summary["outputs"][NOMINAL_TAG] = AnalysisManager.write_events(events, config["variables_of_interest"], output_name, NOMINAL_TAG)
    job_summary["n_events_selected"][NOMINAL_TAG] = len(events)
    t_elapsed = time.time() - t_start

    # Calculate performance metrics
    job_summary["time"] = t_elapsed
    job_summary["time_frac_load"] = t_elapsed_load / t_elapsed
    job_summary["time_frac_samples"] = t_elapsed_samples / t_elapsed
    job_summary["time_frac_syst"] = t_elapsed_syst / t_elapsed
    job_summary["time_frac_taggers"] = t_elapsed_taggers / t_elapsed
    job_summary["successful"] = True
    if not os.path.exists(config["dir"]):
        output_dir = ""
        config["summary_file"] = os.path.split(config["summary_file"])[-1] 
    else: 
        output_dir = os.path.abspath(config["dir"]) + "/"
    output_name = output_dir + config["output_name"]
    # Dump json summary
    if condor:
      with open(config["summary_file"].split("/")[-1], "w") as f_out:
        # json.dump(job_summary, f_out, sort_keys = True, indent = 4)
        f_out.write(orjson.dumps(job_summary, option=orjson.OPT_SERIALIZE_NUMPY).decode("utf-8"))    
    else:
      with open(config["summary_file"], "w") as f_out:
        # json.dump(job_summary, f_out, sort_keys = True, indent = 4)
        f_out.write(orjson.dumps(job_summary, option=orjson.OPT_SERIALIZE_NUMPY).decode("utf-8"))
    return job_summary

class AnalysisManager():
    """
    Manages the running of an entire analysis.
    See HiggsDNA tutorial for more details:
        - https://sam-may.github.io/higgs_dna_tutorial.github.io/

    The ``AnalysisManager`` class is designed such that you can launch an instance through ``scripts/run_analysis.py``, kill the script, and then relaunch to resume running your analysis.
    ``AnalysisManager`` will periodically pickle itself to enable saving/loading its progress.
    This is safe to do and you don't need to worry about jobs being forgotten about and/or being submitted multiple times.
    Most attributes are not allowed to change when loading to prevent things like inconsistent mappings between input files and jobs between runnings. 
    
    If you are running an entirely new analysis, you must specify a ``config`` json or dict with the fields listed in ``REQUIRED_FIELDS``.

    If you are resuming a previously paused/killed instance, it is enough to just pass the ``output_dir`` of your previous run.
    When resuming a previously paused/killed instance (i.e. same ``output_dir``), it is possible to change the attributes listed in ``OVERWRITABLE_FIELDS``.
    All other attributes will be taken from the pickled instance.

    :param output_dir: path to output directory
    :type output_dir: str
    :param config: json/dictionary that specifies the physics content (samples, tag sequence, systematics) of an entire analysis
    :type config: dict or str, optional if and only if you specify an ``output_dir`` that has a pickled ``AnalysisManager`` instance from a previous run
    """
    REQUIRED_FIELDS = ["variables_of_interest", "tag_sequence", "systematics", "samples", "branches"] # these fields **must** be present in config file
    OVERWRITABLE_FIELDS = ["merge_outputs", "unretire_jobs", "retire_jobs", "reconfigure_jobs", "n_cores", "log-level", "log-file", "short", "with_skimmed"] # these fields can be safely changed between different runs
    # TODO: add functionality for switching between batch_system = "local" and batch_system = "condor" between runs
    DEFAULTS = { # kwargs with default values
            "name" : "my_analysis", # doesn't affect actual code, just for more informative printouts
            "function" : {
                "module_name" : "higgs_dna.analysis",
                "function_name" : "run_analysis"
            },
            "batch_system" : "condor",
            #"batch_system" : "local",
            "fpo" : None, # number of input files per output file (i.e. per job)
            "n_cores" : 5, # number of cores for local running
            "merge_outputs" : False,
            "unretire_jobs" : False,
            "retire_jobs" : False,
            "reconfigure_jobs" : False,
            "short" : False, # run only 1 job per sample x year
            "with_skimmed": False # whether to process skimmed files
    }

    def __init__(self, output_dir = "output", config = {}, **kwargs):
        self.output_dir = os.path.abspath(output_dir)
        self.pickle_file = self.output_dir + "/analysis_manager.pkl"
        self.summary_file = self.output_dir + "/summary.json"

        # Resuming a previous run?
        #   Check if there is a pickled AnalysisManager instance. If so, you only need to properly specify ``output_dir`` and everything else will be loaded from the pkl.
        if os.path.exists(self.pickle_file): # resuming a previous run
            logger.warning("[AnalysisManager : __init__] We are loading a saved pickle file '%s' -- please be sure this behavior is inteneded. If you want to run a completely new analysis, please specify a new ``output_dir``." % (self.pickle_file))
            if config:
                logger.warning("[AnalysisManager : __init__] A non-empty ``config`` arg was passed. This will be ignored and the ``AnalysisManager`` state will be loaded from the pkl file. If you are passing the same ``config`` that was used previously, it is safe to ignore this warning.")
            
            # Load pickled ``AnalysisManager`` instance
            with open(self.pickle_file, "rb") as f_in:
                saved_analysis_manager = dill.load(f_in)
            for attr, value in saved_analysis_manager.__dict__.items():
                logger.debug("[AnalysisManager : __init__] Setting attribute '%s' as : %s" % (attr, str(value)))
                setattr(self, attr, value)
            
            # Update overwritable kwargs and warn about any others
            for k,v in kwargs.items():
                if k in self.OVERWRITABLE_FIELDS:
                    # Check if this is already present in pickled instance, and if so, let the user know if it is being changed from its previous value.
                    if hasattr(self, k):
                        if not v == getattr(self, k):
                            logger.info("[AnalysisManager : __init__] Attribute '%s' is being updated from its previous value of '%s' -> '%s'." % (k, getattr(self, k), v))
                    logger.debug("[AnalysisManager : __init__] Setting attribute '%s' as : '%s'." % (k, v)) 
                    setattr(self, k, v)
                elif hasattr(self, k):
                    logger.warning("[AnalysisManager : __init__] kwarg '%s' with value '%s' was given to constructor, but will be ignored in favor of the pickled value '%s'." % (k, v, getattr(self, k)))
                else:
                    logger.warning("[AnalysisManager : __init__] Not sure what to do with kwarg '%s' with value '%s'." % (k, v))

            # Modify ``JobsManager`` and ``Task`` objects based on any updated values for ``OVERWRITABLE_FIELDS``
            self.modify_jobs()


        # Running a new analysis
        else:
            self.config = copy.deepcopy(load_config(config))
            logger.info("[AnalysisManager : __init__] Initializing a new AnalysisManager instance, as no previous pkl file was found to load state from.")

            # Check for required fields in config
            if any([x not in config.keys() for x in self.REQUIRED_FIELDS]):
                logger.exception("[AnalysisManager : __init__] The fields '%s' are required to be present in the ``config`` json/dict, but one or more were not found. The found keys were : '%s'." % (str(self.REQUIRED_FIELDS), str(config.keys())))
                raise ValueError()
            # Set config and kwargs as attributes
            for dicts in [config, kwargs]:
                for k,v in dicts.items():
                    logger.debug("[AnalysisManager : __init__] Setting attribute '%s' as : '%s'." % (k, v))
                    setattr(self, k, v)
            # Set any remaining attributes from default values
            for k,v in self.DEFAULTS.items():
                if not hasattr(self, k):
                    setattr(self, k, v)

            self.update_samples()

            # make output dir
            os.system("mkdir -p %s" % (self.output_dir))

            if self.variables_of_interest is None:
                self.variables_of_interest = []
            
            
            # Create jobs manager
            if self.batch_system.lower() in ["condor", "htcondor"]:
                self.jobs_manager = CondorManager(output_dir = self.output_dir)
            elif self.batch_system.lower() == "local":
                self.jobs_manager = LocalManager(output_dir = self.output_dir, n_cores = self.n_cores)

            # Check if --no_systematics option was given
            if self.no_systematics:
                self.systematics["no_systematics"] = True

            # Create samples manager
            self.sample_manager = SampleManager(**self.samples)

            # Test construction of tag sequence and systematics producer
            self.test_construction()

            self.prepared_analysis = False


    def update_samples(self):
        """
        Propagate correct behavior of --sample_list and --years args from command line.
        """
        # Check for command line updates to samples
        if hasattr(self, "years"):
            years = self.years.split(",")
            if not years == self.samples["years"]:
                logger.warning("[AnalysisManager : update_samples] Years were provided through the config as '%s', but were also specified from the command line as '%s', which is the version we will use." % (str(self.samples["years"]), str(years)))
                self.samples["years"] = years

        if hasattr(self, "sample_list"):
            samples = self.sample_list.split(",")
            if not samples == self.samples["sample_list"]:
                logger.warning("[AnalysisManager : update_samples] Sample list was provided through the config as '%s', but was also specified from the command line as '%s', which is the version we will use." % (str(self.samples["sample_list"]), str(samples)))
                self.samples["sample_list"] = samples

    def modify_jobs(self):
        """
        Modify the behavior of ``JobsManager`` and ``Task`` instances based on new values given through kwargs.
        """
        # Check if we switched batch system
        if self.batch_system == "local" and not isinstance(self.jobs_manager, LocalManager):
            logger.info("[AnalysisManager : modify_jobs] Converting all unfinished jobs from CondorJob -> LocalJob.")
            self.jobs_manager = self.jobs_manager.convert_to_local() # convert from CondorManager -> LocalManager
        elif self.batch_system.lower() in ["condor", "HTCondor"] and not isinstance(self.jobs_manager, CondorManager):
            logger.info("[AnalysisManager : modify_jobs] Converting all unfinished jobs from LocalJob -> CondorJob.")
            self.jobs_manager = self.jobs_manager.convert_to_condor() # convert from LocalManager -> CondorManager
            

        # Update n_cores for ``LocalManager``
        if isinstance(self.jobs_manager, LocalManager):
            self.jobs_manager.n_cores = self.n_cores

        # Check if we were previously running with --short, but was removed for this run
        if not self.short:
            prev_short = any([task.max_jobs >= 0 for task in self.jobs_manager.tasks])
            if prev_short:
                logger.info("[AnalysisManager : modify_jobs] It appears you previously ran with the ``--short`` option but have now removed it. Will submit the full set of jobs for each task.")
                for task in self.jobs_manager.tasks:
                    task.max_jobs = -1
                    task.create_jobs() # create the new jobs (old ones will not be overwritten)
                    self.jobs_manager.customized_jobs = False # need to configure the new jobs
                self.jobs_manager.remerge = self.merge_outputs

        # Reconfigure jobs
        if self.reconfigure_jobs:
            logger.info("[AnalysisManager : modify_jobs] Forcing reconfiguration of jobs (rewriting python configs, executables and condor_submit files).")
            self.jobs_manager.customized_jobs = False
            self.jobs_manager.submit_jobs(summarize = False, dry_run = True)

        # Retire jobs
        if self.retire_jobs:
            logger.warning("[AnalysisManager : modify_jobs] Retiring all unfinished jobs.")
            for task in self.jobs_manager.tasks:
                task.retire_jobs()

        # If ``merge_outputs`` and ``unretire_jobs`` selected, check if any jobs actually get unretired. If so, we need to remerge.
        if self.unretire_jobs:
            logger.info("[AnalysisManager : modify_jobs] Unretiring jobs that failed up to the maximum number of retries.")
            unretired_jobs = False
            for task in self.jobs_manager.tasks:
                task_unretired_jobs = task.unretire_jobs()
                unretired_jobs = unretired_jobs or task_unretired_jobs 

            self.jobs_manager.remerge = unretired_jobs and self.merge_outputs


    def test_construction(self):
        """
        Construct, but do not run, TagSequence and SystematicsProducer.
        """
        logger.info("[AnalysisManager: test_construction] Testing TagSequence and SystematicsProducer construction.")
        test_syst = SystematicsProducer(options = self.systematics)
        test_tags = TagSequence(tag_list = self.tag_sequence)
        logger.info("[AnalysisManager: test_construction] Constructed TagSeqeunce and SystematicsProducer successfully.")


    def run(self):
        """
        Prepare and run analysis until finished, then merge outputs if requested and print out some diagnostic info.
        The AnalysisManager frequently pickles itself throughout this function, so it can easily pick up right where it
        left off if the script is paused or killed.
        """
        logger.debug("[AnalysisManager : run] Running analysis '%s'." % self.name)

        # Measure runtime
        start = time.time()

        # Load samples
        logger.debug("[AnalysisManager : run] Loading samples.")
        self.samples = self.sample_manager.get_samples()
        self.save()

        if not self.prepared_analysis:
            self.prepare_analysis()
            self.save()

        logger.info("[AnalysisManager : run] Running %d tasks with %d total input files split over %d total jobs." % (len(self.jobs_manager.tasks), sum([len(x.files) for x in self.jobs_manager.tasks]), sum([len(x.jobs) for x in self.jobs_manager.tasks])))

        summary = self.jobs_manager.submit_jobs()
        while not self.jobs_manager.complete():
            self.save()
            summary = self.jobs_manager.submit_jobs()

        if self.merge_outputs:
            self.jobs_manager.merge_outputs(self.output_dir)
            self.save()

        self.jobs_manager.summarize()

        end = time.time() - start
        logger.info("[AnalysisManager : run] Finished running analysis '%s'. Elapsed time: %s (hours:minutes:seconds)." % (self.name, str(datetime.timedelta(seconds = end))))

        for task, info in summary.items():
            logger.debug("\n[AnalysisManager : run] Task '%s', PERFORMANCE summary" % (task))
            logger.debug("\t [PERFORMANCE : %s] Processed %d total events in %s (hours:minutes:seconds) of total runtime (%.2f Hz)." % (task, info["physics"]["n_events_initial"], str(datetime.timedelta(seconds = info["performance"]["time"])), float(info["physics"]["n_events_initial"]) / info["performance"]["time"] ))
            for portion in ["load", "syst", "taggers"]:
                if not info["performance"]["time"] > 0:
                    continue
                logger.debug("\t [PERFORMANCE : %s] Fraction of runtime spent on %s: %.2f percent" % (task, portion, (100. * info["performance"]["time_%s" % portion]) / info["performance"]["time"]))

            logger.debug("[AnalysisManager : run] Task '%s', PHYSICS summary" % (task)) 
            for syst_tag, n_events in info["physics"]["n_events_selected"].items():
                if not info["physics"]["n_events_initial"] > 0:
                    continue
                logger.debug("\t [PHYSICS : %s] events set '%s' has eff. of %.2f percent (all taggers)" % (task, syst_tag, (float(n_events) * 100.) / float(info["physics"]["n_events_initial"])))
            if "scale1fb" in info["physics"].keys():
                logger.debug("\t [PHYSICS : %s] With a cross section times BF of %.9f pb and a sum of weights of %.9f, scale1fb for this sample is %.9f" % (task, info["physics"]["norm_factor"], info["physics"]["sum_weights"], info["physics"]["scale1fb"]))


        retired_jobs = []
        for task in self.jobs_manager.tasks:
            retired_jobs += [job for job in task.jobs if job.status == "retired"]

        if retired_jobs:
            logger.warning("[AnalysisManager : run] WARNING: there were %d retired jobs (meaning they failed up to the maximum number of retries). You may want to look into these:" % (len(retired_jobs)))
            for job in retired_jobs:
                logger.warning("\t %s" % job.name_full)

        self.summarize()


    def prepare_analysis(self):
        for sample in self.samples:
            task_branches = self.branches
            # Add extra branches which are always loaded, even if they are not explicitly specified in config
            if sample.is_data:
                task_branches += BRANCHES["data"][sample.year] + BRANCHES["data"]["any"]
            else:
                task_branches += BRANCHES["mc"][sample.year] + BRANCHES["mc"]["any"]

            # Make Sample instance json-able so we can save it in job config files
            jsonable_sample = copy.deepcopy(sample)
            jsonable_sample.files = [x.__dict__ for x in jsonable_sample.files]
            jsonable_sample = jsonable_sample.__dict__
            config = {
                "sample" : copy.deepcopy(jsonable_sample),
                "branches" : task_branches,
                "with_skimmed": self.with_skimmed # Pass the flag to the task config
            }
            # Add skimmed_files to task config if with_skimmed is True
            if self.with_skimmed:
                config["skimmed_files"] = [f.name for f in sample.skimmed_files] if hasattr(sample, "skimmed_files") and sample.skimmed_files else []
            else:
                config["skimmed_files"] = []


            for x in ["systematics", "tag_sequence", "function", "variables_of_interest"]:
                config[x] = copy.deepcopy(getattr(self, x))

            if "branch_map" in self.config.keys():
                config["branch_map"] = copy.deepcopy(self.config["branch_map"])

            self.jobs_manager.add_task(
                    Task(
                        name = sample.name,
                        output_dir = self.output_dir + "/" + sample.name,
                        files = sample.files,
                        config = config,
                        fpo = self.fpo if self.fpo is not None else sample.fpo, 
                        max_jobs = 1 if self.short else -1
                    )
            )
        self.prepared_analysis = True


    def save(self):
        with open(self.pickle_file.replace(".pkl", "_temp.pkl"), "wb") as f_out:
            dill.dump(self, f_out)

        os.system("mv %s %s" % (self.pickle_file.replace(".pkl", "_temp.pkl"), self.pickle_file))


    def summarize(self):
        """
        Save map of sample_name : process_id so different samples can be distinguished in merged outputs.
        TODO: aggregate diagnostic info from Taggers and Systematics 
        """
        self.summary = {
            "sample_id_map" : self.sample_manager.process_id_map,
            "config" : self.config
        }

        for k,v in vars(self).items():
            if k == "summary":
                continue
            if k in self.summary["config"].keys():
                continue # don't double save things
            if is_json_serializable(v):
                self.summary[k] = v

        with open(self.summary_file, "w") as f_out:
            json.dump(self.summary, f_out, sort_keys = True, indent = 4)

    # def get_file_handler(file):
    #     xrootd_src = file.startswith("root://")
    #     if not xrootd_src:
    #         return {"file_handler": uproot.MultithreadedFileSource} # otherwise the memory maps overload available Vmem
    #     elif xrootd_src:
    #         # uncomment below for MultithreadedXRootDSource
    #         return {"xrootd_handler": uproot.source.xrootd.MultithreadedXRootDSource}
    #     # return {}
    #  **get_file_handler(file),

    @staticmethod
    def get_file_handler(filename):
        xrootd_src = filename.startswith("root://")
        if not xrootd_src:
            return {"file_handler": uproot.MultithreadedFileSource} # otherwise the memory maps overload available Vmem
        elif xrootd_src:
            # uncomment below for MultithreadedXRootDSource
            return {"xrootd_handler": uproot.source.xrootd.MultithreadedXRootDSource}
        return {}

    @staticmethod
    def load_events(config):
        """
        Load all branches in ``branches`` from "Events" tree from all nanoAODs in ``files`` into a single zipped ``awkward.Array``.
        Also calculates and returns the sum of weights from nanoAOD "Runs" tree.        

        :param files: list of files
        :type files: list of str
        :param branches: list of branches. If any branches are requested that are not present in the input nanoAOD, they will be silently omitted.
        :type branches: list of str or tuple
        :returns: array of events, sum of weights
        :rtype: awkward.Array, float
        """
        events = []
        sum_weights = 0

        files = config["files"]
        # skimmed_files will be fetched conditionally
        branches = config["branches"]
        is_data = config["sample"]["is_data"]
        year = config["sample"]["year"]

        with_skimmed = config.get("with_skimmed", False)
        events_skimmed_file_for_merging = None # Will hold concatenated events from skimmed files

        if with_skimmed:
            skimmed_files_paths = config.get("skimmed_files", []) # Get paths, default to empty list
            if skimmed_files_paths and not is_data: # Process only if paths exist and not data
                temp_events_skimmed_list = [] # To collect events from each skimmed file
                for skimmed_file in skimmed_files_paths: # Iterate over the provided paths
                    try:
                        f = uproot.open(skimmed_file, timeout = 300, num_workers=1)
                    except Exception:
                        if (os.system(f"xrdcp '{skimmed_file}' '/tmp/{os.getpid()}/{os.path.basename(skimmed_file)}'")):
                            raise RuntimeError("xrdcp failed")
                        f = uproot.open(f'/tmp/{os.getpid()}/{os.path.basename(skimmed_file)}',timeout = 300, num_workers=1)
                    tree = f["Events"]
                    if "Generator_weight" in tree.keys():
                        sum_genWeight = numpy.sum(tree["Generator_weight"])
                        sum_weights += sum_genWeight #FIXME: use genWeight or Generator_weight
                        logger.debug("[AnalysisManager : GeneratorWeightSum] Sum of Generator_weight: {}".format(sum_genWeight))
                        unique_values = numpy.unique(tree["Generator_weight"])
                        for i, value in enumerate(unique_values):
                            unique_counts = numpy.sum(tree["Generator_weight"] == value)
                            logger.debug("[AnalysisManager : GeneratorWeight] Unique values of Generator_weight: {}, numbers: {}".format(value, unique_counts))
                    if is_data:
                        duplicated_sample_remover = DuplicatedSamplesTagger(is_data=is_data)
                        duplicated_remove_cut = duplicated_sample_remover.calculate_selection(skimmed_file, tree)

                        trimmed_branches = [x for x in branches if x in tree.keys()]
                        events_skimmed_file = tree.arrays(trimmed_branches, library = "ak", how = "zip")

                        events_skimmed_file = events_skimmed_file[duplicated_remove_cut]
                    else:
                        mc_overlap_remover = MCOverlapTagger(is_data=is_data, year=year)
                        overlap_cut = mc_overlap_remover.overlap_selection(skimmed_file, tree)

                        trimmed_branches = [x for x in branches if x in tree.keys()]
                        events_skimmed_file = tree.arrays(trimmed_branches, library = "ak", how = "zip") 
                        events_skimmed_file = events_skimmed_file[overlap_cut]
                    f.close()
                    logger.debug("[AnalysisManager : Load events_skimmed_file] Sample type for %s: %s" % (skimmed_file, events_skimmed_file.type))
                    temp_events_skimmed_list.append(events_skimmed_file) # Add processed events to list

                # After processing all skimmed files, if any, concatenate them
                if temp_events_skimmed_list:
                    events_skimmed_file_for_merging = awkward.concatenate(temp_events_skimmed_list)
                    logger.debug("[AnalysisManager : Concatenated skimmed events] Type: %s, Length: %d" % (events_skimmed_file_for_merging.type, len(events_skimmed_file_for_merging)))

            
        for file in files:
            try:
                f = uproot.open(file, timeout = 300, num_workers=1)
            except Exception:
                if (os.system(f"xrdcp '{file}' '/tmp/{os.getpid()}/{os.path.basename(file)}'")):
                    raise RuntimeError("xrdcp failed")
                f = uproot.open(f'/tmp/{os.getpid()}/{os.path.basename(file)}',timeout = 300, num_workers=1)

            runs = f["Runs"]
            # if "genEventCount" in runs.keys() and "genEventSumw" in runs.keys():
            #     # sum_weights += numpy.sum(runs["genEventSumw"].array()) #FIXME: use genWeight or Generator_weight
            #     sum_genWeight = numpy.sum(runs["genEventSumw"].array())
            #     logger.debug("[AnalysisManager : genEventSumw] genEventSumw: {}".format(sum_genWeight))
            # elif "genEventCount_" in runs.keys() and "genEventSumw_" in runs.keys():
            #     # sum_weights += numpy.sum(runs["genEventSumw_"].array()) #FIXME: use genWeight or Generator_weight
            #     sum_genWeight = numpy.sum(runs["genEventSumw_"].array())
            #     logger.debug("[AnalysisManager : genEventSumw_] genEventSumw_: {}".format(sum_genWeight))
            tree = f["Events"]


            if "Generator_weight" in tree.keys():
                sum_genWeight = numpy.sum(tree["Generator_weight"])
                sum_weights += sum_genWeight #FIXME: use genWeight or Generator_weight
                logger.debug("[AnalysisManager : GeneratorWeightSum] Sum of Generator_weight: {}".format(sum_genWeight))
                unique_values = numpy.unique(tree["Generator_weight"])
                for i, value in enumerate(unique_values):
                    unique_counts = numpy.sum(tree["Generator_weight"] == value)
                    logger.debug("[AnalysisManager : GeneratorWeight] Unique values of Generator_weight: {}, numbers: {}".format(value, unique_counts))

            # if "genWeight" in tree.keys():
            #     sum_genWeight = numpy.sum(tree["genWeight"])
            #     # sum_weights += sum_genWeight #FIXME: use genWeight or Generator_weight
            #     logger.debug("[AnalysisManager : GenWeightSum] Sum of genWeight: {}".format(sum_genWeight))
            #     unique_values = numpy.unique(tree["genWeight"])
            #     for value in unique_values:
            #         logger.debug("[AnalysisManager : GenWeight] Unique values of genWeight: ".format(value))
                    
            # Get events that is not duplicated
            if is_data:
                duplicated_sample_remover = DuplicatedSamplesTagger(is_data=is_data)
                duplicated_remove_cut = duplicated_sample_remover.calculate_selection(file, tree)

                trimmed_branches = [x for x in branches if x in tree.keys()]
                # event_file = awkward.Array([])
                # for array in tree.iterate(trimmed_branches, library="ak", how='zip', step_size=100000):
                #     event_file.concatenate(array)
                events_file = tree.arrays(trimmed_branches, library = "ak", how = "zip") #TODO: There is a bug here.

                events_file = events_file[duplicated_remove_cut]
            else:
                mc_overlap_remover = MCOverlapTagger(is_data=is_data, year=year)
                overlap_cut = mc_overlap_remover.overlap_selection(file, tree)

                trimmed_branches = [x for x in branches if x in tree.keys()]
                # event_file = awkward.Array([])
                # for array in tree.iterate(trimmed_branches, library="ak", how='zip', step_size=100000):
                #     event_file.concatenate(array)
                events_file = tree.arrays(trimmed_branches, library = "ak", how = "zip") #TODO: There is a bug here.

                # commented by Pei-Zhu
                events_file = events_file[overlap_cut]

                if int(year[:4]) > 2020:
                    events_file = photon_scale_smear_run3(events_file, year)
                    events_file = electron_scale_smear_run3(events_file, year)

                    # for field in events_file["Photon"].fields:
                    #     logger.info("[AnalysisManager : photon_scale_smear_run3] Photon field: %s" % field)

                if with_skimmed and not is_data and events_skimmed_file_for_merging is not None:

                    # add muon sys
                    events_keys_muon = events_file.Muon.fields  
                    skimmed_keys_muon = events_skimmed_file_for_merging.Muon.fields
                    extra_keys = [key for key in skimmed_keys_muon if key not in events_keys_muon] 
                    for key in extra_keys:
                        events_file['Muon'] = awkward.with_field(events_file['Muon'], events_skimmed_file_for_merging['Muon'][key], key)
                    # add photon sys
                    if (int(year[:4]) < 2020):
                        events_keys_photon = events_file.Photon.fields  
                        skimmed_keys_photon = events_skimmed_file_for_merging.Photon.fields
                        extra_keys = [key for key in skimmed_keys_photon if key not in events_keys_photon] 
                        for key in extra_keys:
                            events_file['Photon'] = awkward.with_field(events_file['Photon'], events_skimmed_file_for_merging['Photon'][key], key)
                    # add jets sys
                    events_keys_jet = events_file.Jet.fields
                    skimmed_keys_jet = events_skimmed_file_for_merging.Jet.fields
                    extra_keys = [key for key in skimmed_keys_jet if key not in events_keys_jet]
                    for key in extra_keys:
                        events_file['Jet'] = awkward.with_field(events_file['Jet'], events_skimmed_file_for_merging['Jet'][key], key)      
                    # add MET sys
                    events_keys_other = events_file.fields
                    skimmed_keys_other = events_skimmed_file_for_merging.fields
                    extra_keys = [key for key in skimmed_keys_other if key not in events_keys_other]
                    for key in extra_keys:
                        events_file = awkward.with_field(events_file, events_skimmed_file_for_merging[key], key)
                    
            f.close()
            
            # # FIXME: DANGEROUS!
            # events_file = events_file[(events_file["run"]==316470) & (events_file["luminosityBlock"]==370) & (events_file["event"]==486186232)]

            logger.debug("[AnalysisManager : Load samples] Sample type: %s" % events_file.type)

            events.append(events_file)

            logger.debug("[AnalysisManager : load_events] Loaded %d events from file '%s'." % (len(events_file), file))

        events = awkward.concatenate(events)

        return events, float(sum_weights)


    @staticmethod
    def write_events(events, save_branches, name, syst_tag):
        """
        For each set of events in ``events``, saves all fields in ``save_branches`` to a .parquet file.

        :param events: dictionary of systematic variations with independent collections : array of events
        :type events: dict
        :param save_branches: list of fields to save in output file
        :type save_branches: list of str or tuple or list
        :param name: name template for output files, which will be updated based on each key of ``events``
        :type name: str
        :returns: dictionary of keys from ``events`` : parquet file
        :rtype: dict
        """
        out_name = "%s_%s.parquet" % (name, syst_tag)

        if not len(events) >= 1:
            return out_name 
        save_map = {}
        for branch in save_branches:
            if isinstance(branch, tuple) or isinstance(branch, list):
                save_name = "_".join(branch)
                if isinstance(branch, list):
                    branch = tuple(branch)
            else:
                save_name = branch
            if isinstance(branch, tuple):
                present = branch[1] in events[branch[0]].fields
            else:
                present = branch in events.fields
            if not present:
                logger.warning("[AnalysisManager : write_events] Branch '%s' was not found in events array. This may be expected (e.g. gen info for a data file), but please ensure this makes sense to you." % str(branch))
                continue
            save_map[save_name] = events[branch]

        for field in events.fields:
            if "weight_" in field and not field in save_map.keys():
                save_map[field] = events[field]

        events = awkward.zip(save_map, depth_limit=1)    
        logger.debug("[AnalysisManager : write_events] Writing output file '%s'." % (out_name))
        if condor:
       	  #awkward.to_parquet(events, out_name.split("/")[-1],list_to32=True) 
            awkward.to_parquet(events, out_name.split("/")[-1]) 
        else:
          #awkward.to_parquet(events, out_name,list_to32=True) 
            awkward.to_parquet(events, out_name) 
        return out_name
