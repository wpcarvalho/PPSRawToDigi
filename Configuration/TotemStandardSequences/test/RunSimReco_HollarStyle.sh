
#!/bin/sh
echo "start submit script"

cmsRun SimStep_HollarStyle.py
cmsRun BTV-Spring14dr-00021_1_cfg.py 
cmsRun BTV-Spring14dr-00021_2_cfg.py
