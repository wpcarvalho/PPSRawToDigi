
#!/bin/sh
echo "start submit script"

cmsRun SimStep_HollarStyle_XML_T1T2CMS.py
cmsRun BTV-Spring14dr-00021_1_cfg.py
cmsRun BTV-Spring14dr-00021_2_cfg.py
cmsRun step2a_CMS_NTUPLE.py