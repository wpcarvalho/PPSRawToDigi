import FWCore.ParameterSet.Config as cms

from EventFilter.CTPPSRawToDigi.diamondVFATRawToDigi_cfi import diamondVFATRawToDigi

diamondRPRawToDigi = diamondVFATRawToDigi.copy()
diamondRPRawToDigi.subSystem = "RP"
