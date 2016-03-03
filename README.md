# totem-offline software

Totem offline software bases on CMSSW framework and provides plugins that 
simulate and reconstruct beam interactions in Totem experiment. 


## release/7.5.0

Release 7.5.0 contains some of modules from release/7.0.4, responsible for simulation, 
migrated to newer version of cms framework.

### Getting started

```
ssh -X $USER@lxplus
export SCRAM_ARCH=slc6_amd64_gcc491
source /afs/cern.ch/cms/cmsset_default.sh
scram project CMSSW CMSSW_7_5_0
git clone ssh://git@gitlab.cern.ch:7999/totem/totem-offline.git CMSSW_7_5_0/src
cd CMSSW_7_5_0/src
git checkout release/7.5.0
cmsenv
scram build -j 15
```

### Running sample configuration

PPS simualtion config:
```
cmsRun src/test.py
```
