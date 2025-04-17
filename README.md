This is a tool for analyzing AOD files and producing L1T muon ntuples.  

**NOTE:** Before setting up, check the CMSSW release used to create the AOD file you want to process. Older releases may not work! You can find the release in [DAS](https://cmsweb.cern.ch/das/), e.g.: 
```
release dataset=/Cosmics/Run2024I-PromptReco-v1/AOD
``

Setup
-----
```
cmsrel CMSSW_14_0_16
cd CMSSW_14_0_16/src
cmsenv
voms-proxy-init --voms cms --valid 24:00:00

https://github.com/muonDPG/MuonAODAnalyzer.git
scram b -j8
```

Run
-----
To run locally on 1 file (testing):
```
cmsRun L1TMuonAODAnalyzer/pset_files/run3_data.py
```
To submit to CRAB:
```
cd L1TMuonAODAnalyzer/pset_files/
crab submit
```

Collections
-----------
- ```l1t::Muon``` objects  
- ```reco::Muon``` objects  
- ```reco::Track``` objects  
- BMTF ```l1t::RegionalMuonCand``` objects   

Useful links
-------------
CRAB tutorial: https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookCRAB3Tutorial  
CRAB config file: https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3ConfigurationFile  
CRAB exit codes: https://twiki.cern.ch/twiki/bin/view/CMSPublic/JobExitCodes  