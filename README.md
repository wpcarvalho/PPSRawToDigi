# totem-offline software

Totem offline software bases on CMSSW framework and provides plugins that 
simulate and reconstruct beam interactions in Totem experiment. 


## release/8.0.X

Release 8.0.X contains some of modules from release/7.0.4, responsible for recunstruction, 
migrated to newer version of cms framework.

### Getting started

```
ssh -X $USER@lxplus
export SCRAM_ARCH=slc6_amd64_gcc493
source /afs/cern.ch/cms/cmsset_default.sh
cmsrel CMSSW_8_0_0_pre5
git clone ssh://git@gitlab.cern.ch:7999/totem/totem-offline.git CMSSW_8_0_0_pre5/src
cd CMSSW_8_0_0_pre5/src
git checkout release/8.0.X
cmsenv
scram build -j 15
```

### Running sample configuration

RP reco RunII config:
```
cmsRun src/config_9998.152127.py
```
