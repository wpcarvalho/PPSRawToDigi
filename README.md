[![Build Status](https://travis-ci.org/CTPPS/ctpps-offline.svg?branch=develop)](https://travis-ci.org/CTPPS/ctpps-offline)

# ctpps-offline software

CTPPS offline software bases on CMSSW framework. 


## release/8.1.X

Release 8.1.X contains some of modules responsible for reconstruction,
migrated to newer version of CMSSW framework.

### Getting started

```
ssh -X $USER@lxplus
source /afs/cern.ch/cms/cmsset_default.sh
cmsrel CMSSW_8_1_0_pre5
git clone https://github.com/setesami/PPSRawToDigi.git CMSSW_8_1_0_pre5/src
cd CMSSW_8_1_0_pre5/src
cmsenv
scram build -j 15
```

### Running sample configuration

```
cmsRun src/raw_data_chain_test.py

### Runing samples config for diamond

cmsRun src/diamondraw_data_chain_test.py
```
