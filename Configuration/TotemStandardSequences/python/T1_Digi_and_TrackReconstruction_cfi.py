# No pile up for the mixing module
#from SimGeneral.MixingModule.mixNoPU_cfi  import *
from Configuration.TotemCommon.mixNoPU_cfi import *
########################### DIGI + RECO T1 ##########################################

from SimTotem.T1Digitizer.T1DigisVFAT_cfi import *

from RecoTotemT1T2.T1MakeCluster.T1MakeCluster_cfi import *

from RecoTotemT1T2.T1RecHit.T1RecHit_cfi import *

from RecoTotemT1T2.T1RoadProducer.T1RoadProducer_cfi import *

from RecoTotemT1T2.T1TrackProducer2.T1TrackProducer2_cfi import *
