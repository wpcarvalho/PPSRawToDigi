import FWCore.ParameterSet.Config as cms

XMLIdealGeometryESSource = cms.ESSource("XMLIdealGeometryESSource",
    geomXMLFiles = cms.vstring('Geometry/CMSCommonData/data/materials.xml',
                               'Geometry/CMSCommonData/data/rotations.xml',
                               'Geometry/CMSCommonData/data/extend/cmsextent.xml',
                               'Geometry/CMSCommonData/data/cms.xml',
                               'Geometry/CMSCommonData/data/cmsMother.xml',
                              # 'Geometry/XMLTutorial/data/main.xml'
                               'Geometry/XMLTutorial/data/20160823/CTPPS_Diamond_Materials.xml',
                               'Geometry/XMLTutorial/data/20160823/CTPPS_Diamond_Transformations.xml',
                               'Geometry/XMLTutorial/data/20160823/CTPPS_Diamond_X_Distance.xml',
                               'Geometry/XMLTutorial/data/20160823/CTPPS_Diamond_Parameters.xml',
                               'Geometry/XMLTutorial/data/20160823/DiSeg/Pat1_Str1.xml',
                               'Geometry/XMLTutorial/data/20160823/DiSeg/Pat2_Str1.xml',
                               'Geometry/XMLTutorial/data/20160823/DiSeg/Pat2_Str2.xml',
                               'Geometry/XMLTutorial/data/20160823/DiSeg/Pat3_Str1.xml',
                               'Geometry/XMLTutorial/data/20160823/DiSeg/Pat3_Str2.xml',
                               'Geometry/XMLTutorial/data/20160823/DiSeg/Pat3_Str3.xml',
                               'Geometry/XMLTutorial/data/20160823/DiSeg/Pat3_Str4.xml',
                               'Geometry/XMLTutorial/data/20160823/DiSeg/Pat4_Str1.xml',
                               'Geometry/XMLTutorial/data/20160823/DiSeg/Pat4_Str2.xml',
                               'Geometry/XMLTutorial/data/20160823/DiSeg/Pat4_Str3.xml',
                               'Geometry/XMLTutorial/data/20160823/DiSeg/Pat4_Str4.xml',
                               'Geometry/XMLTutorial/data/20160823/DiSeg/Pat4_Str5.xml',
                               'Geometry/XMLTutorial/data/20160823/DiStAsmb/Plane1.xml',
                               'Geometry/XMLTutorial/data/20160823/DiStAsmb/Plane2.xml',
                               'Geometry/XMLTutorial/data/20160823/DiStAsmb/Plane3.xml',
                               'Geometry/XMLTutorial/data/20160823/DiStAsmb/Plane4.xml',
                               'Geometry/XMLTutorial/data/20160823/CTPPS_Diamond_Detector_Assembly.xml',
                               'Geometry/XMLTutorial/data/20160823/CTPPS_Timing_Horizontal_Pot.xml',
                               'Geometry/XMLTutorial/data/20160823/CTPPS_Timing_Positive_Station.xml',
                               'Geometry/XMLTutorial/data/20160823/CTPPS_Timing_Negative_Station.xml',
                               'Geometry/XMLTutorial/data/20160823/CTPPS_Timing_Stations_Assembly.xml',


),
    rootNodeName = cms.string('cms:OCMS')

)
