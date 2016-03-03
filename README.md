# totem-offline software

Totem offline software bases on CMSSW framework and provides plugins that 
simulate and reconstruct beam interactions in Totem experiment. 


## release/7.0.4

Release 7.0.4 is basically the ctpps branch taken from Totem svn (https://svnweb.cern.ch/cern/wsvn/totem/branches/CMSSW_7_0_4_ctpps/)
with revision 10969. With applied additional changes from branch CMSSW_7_0_4 up to revision 11589.

### Getting started

```
ssh -X $USER@lxplus
export SCRAM_ARCH=slc6_amd64_gcc481
source /afs/cern.ch/cms/cmsset_default.sh
scram project CMSSW CMSSW_7_0_4
git clone ssh://git@gitlab.cern.ch:7999/totem/totem-offline.git CMSSW_7_0_4/src
cd CMSSW_7_0_4/src
git checkout release/7.0.4
cmsenv
scram build -j 15
```

### Running sample configuration

TOTEM RP config:
```
cmsRun src/Configuration/TotemStandardSequences/test/RP/prodRPinelasticBeta90Energy6500GeV_cfg.py
```

TOTEM T1T2 minimum bias config:
```
cmsRun src/Configuration/TotemStandardSequences/test/T1T2/prodT1T2phojetMBbeta90energy6500eV_cfg.py
```

TOTEM RPT1T2 single diffraction config:
```
cmsRun src/Configuration/TotemStandardSequences/test/RPT1T2/prodRPT1T2pythiaSDbeta90energy6500GeV_cfg.py
```
