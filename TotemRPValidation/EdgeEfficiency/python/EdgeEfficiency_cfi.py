import FWCore.ParameterSet.Config as cms

EdgeEfficiency = cms.EDAnalyzer("EdgeEfficiency",
    	Verbosity = cms.uint32(0),
	RootFileName = cms.string('EdgeEfficiency.root'),
	selectedTest = cms.uint32(0),			# LHC run -> 0, RP test beam -> 1, detector package -> 2

	#tolerance_in_v = cms.double(1),		# [mm] a track is accepted if it falls in the radius of a hit with tolerance_in_v
	tolerance_in_v = cms.double(3 * 0.066),		# [mm] a track is accepted if it falls in the radius of a hit with tolerance_in_v
	tolerance_in_angle = cms.double(1E-2), 		# [rad] a track is dropped if the angle is larger than this
	RP_angle_resolution = cms.double(0.000658073),  # [rad]
	RP_size_along_the_z_axis = cms.double(40.6),     # [mm]

	RPFittedTrackCollectionLabel = cms.InputTag("RPSingleTrackCandCollFit"),
	DetSetVectorRPRecoHitLabel = cms.InputTag("RPHecoHitProd")
)
