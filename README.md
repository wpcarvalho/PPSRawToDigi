[![Build Status](https://travis-ci.org/CTPPS/ctpps-offline.svg?branch=develop_80x)](https://travis-ci.org/CTPPS/ctpps-offline)

# ctpps-offline software

CTPPS offline software based on CMSSW framework.


## release/8.0.X

Release 8.0.X contains some of modules responsible for reconstruction,
migrated to newer version of CMSSW framework.

It is backport from release 8.1.X

### Getting started

```
ssh -X $USER@lxplus
source /afs/cern.ch/cms/cmsset_default.sh
cmsrel CMSSW_8_0_8
git clone https://github.com/CTPPS/ctpps-offline.git CMSSW_8_0_8/rc
cd CMSSW_8_0_8/src
cmsenv
scram build -j 15
```

### Running sample configuration

```
cmsRun src/raw_data_chain_test.py
```
