"""
Much of the code is taken from:
    ProjectMetis metis/Task.py: https://github.com/aminnj/ProjectMetis/blob/f9e71556cb84496731fa71dcab2dfc82b6e3022f/metis/Task.py
    Author: Nick Amin
"""
import copy
import os
import math
import json
import awkward
import numpy
from tqdm import tqdm

import logging
# logger = logging.getLogger(__name__)
from higgs_dna.utils.logger_utils import simple_logger
logger = simple_logger(__name__)

from higgs_dna.job_management.jobs import Job 
from higgs_dna.utils import awkward_utils
from higgs_dna.utils.misc_utils import create_chunks
from higgs_dna.utils.progress_bar import ProgressBar
from higgs_dna.constants import CENTRAL_WEIGHT

class Task():
    """
    A ``Task`` splits input files for a given sample x year among jobs and performs the monitoring, bookkeeping, and post-processing. 
    It shares management of ``Job`` activity with the ``JobsManager`` class, which configures the jobs for local/condor submission and handles the actual submission.

    :param name: name to identify this Task
    :type name: str
    :param output_dir: base directory for this Task. Each ``Job`` will have its own subdirectory within.
    :type output_dir: str
    :param files: list of input files for this sample x year
    :type files: list
    :param config: dictionary that specifies the physics content (samples, tag sequence, systematics) of this Task and its subsequent ``Job``s
    :type config: dict
    """
    def __init__(self, name, output_dir, files, config, **kwargs):
        self.name = name
        self.output_dir = os.path.abspath(output_dir)
        self.files = files
        self.config = config
        self.job_config = copy.deepcopy(config)
        self.job_config["sample"]["files"] = None # lightweight config to send to jobs

        self.min_completion_frac = kwargs.get("min_completion_frac", 1.0)
        self.fpo = kwargs.get("fpo", None) # number of input files per job
        self.max_jobs = kwargs.get("max_jobs", None) 

        for key, val in kwargs.items():
            setattr(self, key, val)

        if self.max_jobs is None:
            self.max_jobs = -1

        if self.fpo is None:
            if hasattr(config["sample"], "fpo"):
                self.fpo = config["sample"].fpo
            elif isinstance(config["sample"], dict):
                if "fpo" in config["sample"].keys():
                    self.fpo = config["sample"]["fpo"]
        if self.fpo is None:
            self.fpo = 1
  

        os.system("mkdir -p %s" % self.output_dir)

        self.complete = False
        self.merged_output_files = False
        self.wrote_process_ids = False
        self.wrote_years = False
        self.pbar = ProgressBar(self.name)

        self.phys_summary = {
                "n_events_initial" : 0,
                "n_events_selected" : {},
                "sum_weights" : 0.0,
        }

        self.performance = {
                "time" : 0.0,
                "time_load" : 0.0,
                "time_syst" : 0.0,
                "time_taggers" : 0.0
        }
 

        self.jobs = []
        self.create_jobs()

        self.summary = {}
        self.outputs = {}
        self.merged_outputs = {}
        self.summarize()
        

    def create_jobs(self):
        """
        Split input files across jobs per the specified number of files per job (``fpo``) and create virtual ``Job`` instances.
        """
        logger.info("[Task: create_jobs] Task '%s' : splitting %d input files into %d jobs" % (self.name, len(self.files), math.ceil(float(len(self.files)) / float(self.fpo)))) 
        file_splits = create_chunks(self.files, self.fpo)
        file_splits = file_splits
        for idx, file_split in enumerate(tqdm(file_splits)):
            if self.max_jobs >= 0:
                if idx >= self.max_jobs:
                    continue

            # Check if this job was already previously added (may happen if user runs with --short and then removes it)
            if idx < len(self.jobs):
                continue

            self.jobs.append(
                    Job(
                        name = self.name,
                        function = copy.deepcopy(self.config["function"]),
                        inputs = file_split,
                        dir = self.output_dir,
                        idx = idx + 1,
                        config = copy.deepcopy(self.job_config)
                    )
            )

            self.complete = False # if we just added a new job, not complete yet (had to add this for when running once with --short and then rerunning without it)
            self.merged_output_files = False


    def process(self, job_map = None):
        """
        Monitor job progress and perform post-processing when complete.
        
        :param job_map: a dictionary mapping jobs to their process/cluster ids, provided by the JobsManager
        :type job_map: dict
        """
        if self.complete:
            for syst_tag, output in self.merged_outputs.items():
                if not os.path.exists(output):
                    logger.warning("[Task : process] A file may have been deleted. Task '%s' was marked as complete, but output '%s' is not present." % (self.name, output))
            if not self.merged_output_files:
                self.merge_outputs()
            return

        # Check status of all jobs
        for job in self.jobs:
            if job_map is not None:
                if os.path.exists(job.summary_file):
                    job.status = "completed"
                elif job.cluster_id in job_map.keys():
                    if job_map[job.cluster_id] is not None:
                        job.status = job_map[job.cluster_id]
                        continue
                elif job.status == "retired":
                    continue
                else:
                    job.status = "failed"
            else:
                job.query_status()

        # Are we done?
        self.n_completed_jobs = len([job for job in self.jobs if job.status == "completed"]) 
        self.n_retired_jobs = len([job for job in self.jobs if job.status == "retired"]) 
        self.completion_frac = float(self.n_completed_jobs) / float(len(self.jobs))
        self.completion_or_retired_frac = float(self.n_completed_jobs + self.n_retired_jobs) / float(len(self.jobs))

        # Did we successfully process minimum fraction of completed jobs?
        if self.completion_frac >= self.min_completion_frac:
            logger.info("[Task : process] Task '%s' COMPLETED : %d/%d (%.2f percent) of jobs completed, which is >= the minimum job completion fraction for this task (%.2f percent)." % (self.name, self.n_completed_jobs, len(self.jobs), 100. * self.completion_frac, 100. * self.min_completion_frac))
            self.complete = True # Done

        # Otherwise, are the uncompleted jobs "retired" (meaning they failed up to the maximum number of retries)?
        # If so, we will mark this as done but give you a warning about the retired jobs.
        # You can resubmit these jobs with the option --unretire_jobs in run_analysis.py
        elif self.completion_or_retired_frac >= self.min_completion_frac:
            logger.info("[Task : process] Task '%s' COMPLETED : %d/%d (%.2f percent) of jobs completed/retired which is >= the minimum job completion fraction for this task (%.2f percent)." % (self.name, self.n_completed_jobs + self.n_retired_jobs, len(self.jobs), 100. * self.completion_or_retired_frac, 100. * self.min_completion_frac))
            retired_jobs = [job for job in self.jobs if job.status == "retired"]
            for job in retired_jobs:
                logger.warning("[Task : process] WARNING: Task '%s' had to retire job '%s' since it ran unsuccessfully for %d straight times. If this is an MC sample, this will just reduce your statistics. If this is a data job, you have processed less events than you intended!" % (self.name, job.name_full, job.n_attempts))
            self.complete = True

        # If neither of the first three, we are not done yet
        else:
            self.summarize()
            return
        logger.info("[Task : process] Task '%s' : %d/%d (%.2f percent) of jobs completed/retired, which is >= the minimum job completion fraction for this task (%.2f percent). Marking task as complete." % (self.name, self.n_completed_jobs + self.n_retired_jobs, len(self.jobs), 100. * self.completion_or_retired_frac, 100. * self.min_completion_frac))
        # Clean up: kill any idle or running jobs
        n_killed = 0
        jobs_to_kill = [job for job in self.jobs if (job.status == "idle" or job.status == "running")]
        for job in jobs_to_kill:
            logger.debug("[Task : process] Task '%s' : since Task is COMPLETED, we are killing job '%s' with status '%s'" % (self.name, job.name_full, job.status))
            job.kill()
            job.status = "retired"
        jobs_to_retire = [job for job in self.jobs if not job.status == "completed"]
        for job in jobs_to_retire:
            job.force_retirement = True
        if n_killed > 0:
            logger.info("[Task : process] Task '%s' : since Task is COMPLETED, we killed all idle and running jobs (%d jobs killed)" % (self.name, n_killed))
        self.complete = True
        self.summarize()
        return

    
    def unretire_jobs(self):
        """
        Jobs that fail 5 times in a row are permanently "retired" by HiggsDNA, meaning they will not be submitted again.
        For MC, retired jobs just amount to a loss in statistics, as the scale1fb will be calculated based on the sum of event weights for successful jobs.
        For Data, however, there is no such trick and the user may want to submit some data jobs more than 5 times.
        The ``--unretire_jobs`` flag can be given to ``scripts/run_analysis.py`` in order to unretire and resubmit all retired jobs.
        """
        retired_jobs = [job for job in self.jobs if job.status == "retired"]
        if len(retired_jobs) == 0:
            return False

        logger.info("[Task : unretire_jobs] Task '%s' : since you are re-running with option --unretire_jobs, we are resubmitting %d jobs which were previously retired (meaning they repeatedly failed up to the number of max attempts)." % (self.name, len(retired_jobs)))
        for job in retired_jobs:
            logger.debug("[Task : unretire_jobs] Task '%s' : resubmitting job '%s' with status '%s' and %d previously failed attempts." % (self.name, job.name_full, job.status, job.n_attempts))
            job.n_attempts = 0
            job.force_retirement = False
            job.status = "waiting"
    
        self.complete = False
        self.merged_output_files = False
        self.wrote_years = False
        self.wrote_process_ids = False
        self.min_completion_frac = 1.0
        return True

    
    def retire_jobs(self):
        """
        All jobs that are not already complete will be retired.
        This can be undone through the ``unretire_jobs`` method.
        """
        self.process() 
        if self.n_completed_jobs > 0:
            logger.info("[Task : process] Task '%s' : `--retire_jobs` option was selected. %d/%d (%.2f percent) of jobs completed for this task and all others will be retired." % (self.name, self.n_completed_jobs, len(self.jobs), 100. * self.completion_frac)) 
            for job in self.jobs:
                job.force_retirement = True
            self.min_completion_frac = 0.000001

        else:
            logger.info("[Task : process] Task '%s' : `--retire_jobs` option was selected, but no jobs have completed yet for this task. Will set the minimum completion fraction to 0.000001, which will effectively retire all jobs once at least 1 finishes." % (self.name))
            self.min_completion_frac = 0.000001


    def summarize(self):
        """

        """
        # Summarize jobs
        job_summary = {}
        job_summary["all"] = len(self.jobs)
        for status in ["waiting", "idle", "running", "failed", "completed", "retired"]:
            job_summary[status] = len([job for job in self.jobs if job.status == status])
        self.summary["jobs"] = job_summary

        #if not self.complete:
        #    return

        completed_jobs = [job for job in self.jobs if job.status == "completed"]
        for job in completed_jobs:
            if job.processed: # don't record metadata twice
                continue

            # Try to open json summary file for the job
            # Wrapped in a try-except to avoid errors where the file exists but has not finished copying from the remote job
            copied = False
            while not copied:
                try:
                    with open(job.summary_file, "r") as f_in:
                        job_info = json.load(f_in)
                        copied = True
                except:
                    os.system("sleep 0.1s")
            self.phys_summary["n_events_initial"] += job_info["n_events"]
            if not self.config["sample"]["is_data"]:
                self.phys_summary["sum_weights"] += job_info["sum_weights"]

            for syst_tag, n_events in job_info["n_events_selected"].items():
                if syst_tag not in self.phys_summary["n_events_selected"].keys():
                    self.phys_summary["n_events_selected"][syst_tag] = n_events
                else:
                    self.phys_summary["n_events_selected"][syst_tag] += n_events

            for syst_tag, output in job_info["outputs"].items():
                if syst_tag not in self.outputs.keys():
                    self.outputs[syst_tag] = []
                if not os.path.exists(output):
                    if os.path.exists(job.output_dir + "/" + output):
                        output = job.output_dir + "/" + output
                    else:
                        logger.exception("[Task : summarize] Did not find output for job '%s' with dir '%s', output dir '%s', config file '%s', and summary file '%s'." % (job.name, job.dir, job.output_dir, job.config_file, job.summary_file))
                        # raise RuntimeError() # FIXME : this is a temporary fix to avoid crashing, but we should handle this more gracefully
                if not job_info["n_events_selected"][syst_tag] > 0: # skip empty parquet files to avoid errors
                    continue
                else:
                    self.outputs[syst_tag].append(output)

            self.performance["time"] += job_info["time"]
            for portion in ["load", "syst", "taggers"]:
                self.performance["time_%s" % portion] += float(job_info["time"]) * job_info["time_frac_%s" % portion]

            job.processed = True # make sure we don't record its metadata again


        if not self.config["sample"]["is_data"]:
            self.phys_summary["norm_factor"] = self.config["sample"]["norm_factor"]
            if self.phys_summary["sum_weights"] > 0.:
                self.phys_summary["scale1fb"] = (self.config["sample"]["norm_factor"] * 1000.) / self.phys_summary["sum_weights"]
            else:
                self.phys_summary["scale1fb"] = 0.
            self.lumi = self.config["sample"]["lumi"]
            self.scale1fb = self.phys_summary["scale1fb"]
        self.summary["physics"] = self.phys_summary
        self.summary["performance"] = self.performance
        self.pbar.update(job_summary, self.performance, self.phys_summary)

        if self.complete:
            self.merge_outputs()


    def merge_outputs(self):
        """
        Merge output files from all completed jobs into a single file per systematic with independent collection.
        If MC, we also scale the central weight by scale1fb and luminosity, where scale1fb is calculated
        from the sum of weights of completed jobs. We also add a branch CENTRAL_WEIGHT + "_no_lumi" which has just
        scale1fb applied and is not scaled by luminosity.
        """
        self.merged_outputs = {}
        for syst_tag, outputs in self.outputs.items():
            if not outputs:
                continue
            
            merged_output = self.output_dir + "/merged_%s.parquet" % (syst_tag) 
            self.merged_outputs[syst_tag] = merged_output
            merged_events = []

            # FIXME : merging could be improved so that we avoid merging huge numbers of events into a single file and instead split them across multiple files
            for output in outputs:
                merged_events.append(awkward.from_parquet(output))

            logger.debug("[Task : merge_outputs] Task '%s' : merging %d outputs into file '%s'." % (self.name, len(outputs), merged_output))

            merged_events = awkward.concatenate(merged_events)
            if not self.config["sample"]["is_data"]:
                logger.debug("[Task : merge_outputs] Task '%s' : Applying scale1fb and lumi. Scaling central weight branch '%s' in output file '%s' by scale1fb (%.13f) times lumi (%.2f). Adding branch '%s' in output file which has no lumi scaling applied." % (self.name, CENTRAL_WEIGHT, merged_output, self.scale1fb, self.lumi, CENTRAL_WEIGHT + "_no_lumi"))
                #logger.debug(f"[Task : sumWeight is {self.phys_summary['sum_weights']}]")
                #logger.debug(f"[Task : typical central weight is {merged_events['weight_central']}]")

                central_weight = merged_events[CENTRAL_WEIGHT].to_numpy().astype('float64') * self.scale1fb * self.lumi
                central_weight_no_lumi = merged_events[CENTRAL_WEIGHT].to_numpy().astype('float64') * self.scale1fb

                awkward_utils.add_field(
                        events = merged_events,
                        name = CENTRAL_WEIGHT,
                        data = central_weight,
                        overwrite = True
                )
                awkward_utils.add_field(
                        events = merged_events,
                        name = CENTRAL_WEIGHT + "_no_lumi",
                        data = central_weight_no_lumi,
                        overwrite = False # merging and scale1fb application should always be done from unmerged files, which should not have this branch already. If they somehow do have the branch, that is a bad sign something is going wrong...
                )

            awkward.to_parquet(merged_events, merged_output)

        self.wrote_process_ids = False
        self.wrote_years = False
        self.merged_output_files = True


    def add_process_id(self):
        """
        Add a field "process_id" with the process_id value for this sample given by the SampleManager.
        Note that although there are separate ``Sample`` and ``Task`` instances for different years of the same sample,
        the "process_id" field is harmonized across years.
        """
        if self.wrote_process_ids:
            return

        self.process_id = self.config["sample"]["process_id"]
        if self.process_id is None:
            return

        for syst_tag, merged_output in self.merged_outputs.items():
            events = awkward.from_parquet(merged_output)

            if "process_id" in events.fields:
                return

            logger.debug("[Task : add_process_id] Task '%s' : adding field 'process_id' with value %d in output file '%s'" % (self.name, self.process_id, merged_output))
            awkward_utils.add_field(
                    events = events,
                    name = "process_id",
                    data = numpy.ones(len(events)) * self.process_id
            )
            awkward.to_parquet(events, merged_output)
        

        self.wrote_process_ids = True

    
    def add_year(self):
        """
        Add a field "year" with the year for this sample, as taken from the SampleManager.
        Note that because the year is taken from the SampleManager, it is possible the "year" field
        is not the same as the year for the CMS sample (e.g. this might happen on purpose if you were missing a sample 
        for 2016 and used a sample from 2017 in its place).
        The "year" field is stored as an integer.
        """
        if self.wrote_years:
            return

        self.year = self.config["sample"]["year"]

        for syst_tag, merged_output in self.merged_outputs.items():
            events = awkward.from_parquet(merged_output)

            if "year" in events.fields:
                return

            logger.debug("[Task : add_process_id] Task '%s' : adding field 'year' with value '%s' in output file '%s'." % (self.name, self.year, merged_output))

            year_array = numpy.empty(len(events), dtype="S10")
            year_array[:] = self.year
            awkward_utils.add_field(
                    events = events,
                    name = "year",
                    data = year_array 
            )

            awkward.to_parquet(events, merged_output)


        self.wrote_years = True
