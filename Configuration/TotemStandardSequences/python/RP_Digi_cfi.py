import FWCore.ParameterSet.Config as cms

# Use random number generator service
from Configuration.TotemCommon.RandomNumbers_cfi import *

# Use particle table
from SimGeneral.HepPDTESSource.pdt_cfi import *

# No pile up for the mixing module
from SimGeneral.MixingModule.mixNoPU_cfi import *

# RP Strip digitization
from SimTotem.RPDigiProducer.RPSiDetConf_cfi import *
#RPSiDetDigitizer.RPVerbosity = 1