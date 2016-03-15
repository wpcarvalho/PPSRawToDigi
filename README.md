# ctpps-offline software

CTPPPS offline software bases on CMSSW framework. 


## release/8.0.X

Release 8.0.X contains some of modules responsible for recunstruction, 
migrated to newer version of cms framework.

### Getting started

```
ssh -X $USER@lxplus
export SCRAM_ARCH=slc6_amd64_gcc493
source /afs/cern.ch/cms/cmsset_default.sh
cmsrel CMSSW_8_0_0_pre5
git clone https://github.com/CTPPS/ctpps-offline.git CMSSW_8_0_0_pre5/src
cd CMSSW_8_0_0_pre5/src
cmsenv
scram build -j 15
```

### Running sample configuration

RP reco RunII config:
```
cmsRun src/config_9998.152127.py
```
