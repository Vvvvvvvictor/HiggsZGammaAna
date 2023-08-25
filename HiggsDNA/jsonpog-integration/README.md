# jsonPOG-integration
 

**POG folder in gitlab**

In this folder we store all the corrections.
Each object have a json, and each POG has a folder for storage

| directory  | year campaign | name.json |
| ---------- | --------------| ----------|
| POG/EGM  | "X_Y"  | photon.json |
|          |  "" | electron.json |
| POG/TAU  |  "" | tau.json |
| POG/MUON |  "" | muon.json |
| POG/JME  |  "" | fatJetPuppi.json |
|          |  "" | jetCHS.json |
| POG/BTV  |  "" | bjets.json |
|          |  "" | cjets.json |
| POG/LUM  |  "" | puWeights.json | 

Initial notes: 
1. different campaings are organized in folders with label "X_Y" i.e. (2016preVFP_UL, 2016postVFP_UL, 2017_UL, 2018_UL ...) (later we will add the *_EOY)
2. each object in nano get one json
3. the "inputs" labels should be unique and standardized
- (the X_Y year) should be kept in the json to mantain the provenance
    X: 2016preVFP, 2016postVFP, 2017, 2018
    Y: Prompt,EOY, UL
    i.e. 2016postVFP_UL
- second set use the string "sf","syst" (all lower cases)
4. store the json in .gz format for compression


**Repository with templates and tools in github**

for the time being things are here
https://github.com/cms-nanoAOD/correctionlib


**Script folder in gitlab**

Once a PR is made, a test is being started.
The tests will happen with the script defined here
Goal of the test:
* verify can be loaded
* verify that is compliant to the format "Schema Version" 
* run the summary to inspect them 
* at each change we create a verisioned object file and a tag for the whole set Vxx_yy.

**Distribution to users**


_/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration_
is available and synched daily with latest commit of master branch of https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration


**HOW to read the jsons**

Format of the evaluation i.e. for electron
`evaluator[ JSON Map Name ].evaluate( year , Value Type , Working Point , eta , pt )`

See folder Examples on how to read a scale factor value


**INSPECT the jsons**

locally `correction summary file.json`
look what is available here https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/

