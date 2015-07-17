import FWCore.ParameterSet.Config as cms

# beam/optics parameters
from Configuration.TotemOpticsConfiguration.OpticsConfig_4000GeV_0p6_143urad_cfi import *

# geometry
from Configuration.TotemCommon.geometryRP_real_cfi import *
XMLIdealGeometryESSource.geomXMLFiles.append('Geometry/TotemRPData/data/RP_Garage/RP_Dist_Beam_Cent.xml')

from TotemAlignment.RPDataFormats.TotemRPIncludeAlignments_cfi import *

import os

working_dir = os.path.dirname(os.path.realpath(__file__))

alignment_files = cms.vstring(
  working_dir+"/../alignment_station_200m.xml",                # shifts 147m station to 200m, puts all RP edges to the beam
  working_dir+"/rp_position.xml",                              # reasonable RP distance from beam
  working_dir+"/../misalignments/geometry_real_misaligned.xml"      # misalignement
)

TotemRPIncludeAlignments.RealFiles = alignment_files
TotemRPIncludeAlignments.MisalignedFiles = alignment_files

# misalignement measurement error
TotemRPIncludeAlignments.MisalignedFiles.append(working_dir+"/../misalignments/geometry_misaligned.xml")
