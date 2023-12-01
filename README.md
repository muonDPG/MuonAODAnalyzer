This is a tool for analyzing AOD files and producing L1T muon ntuples.  

Setup
-----
```
cmsrel CMSSW_13_0_3
cd CMSSW_13_0_3/src
cmsenv
voms-proxy-init --voms cms --valid 24:00:00

git clone https://github.com/yiannispar/L1TMuonAODAnalyzer.git
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
- ```gmtStage2Digis``` (```l1t::Muon``` objects)  
- ```muons``` (```reco::Muon``` objects)  
- ```standAloneMuons``` (```reco::Track``` objects)  

Useful links
-------------
CRAB tutorial: https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookCRAB3Tutorial  
CRAB config file: https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3ConfigurationFile  
CRAB exit codes: https://twiki.cern.ch/twiki/bin/view/CMSPublic/JobExitCodes  