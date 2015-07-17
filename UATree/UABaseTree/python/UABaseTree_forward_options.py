
# General switches --------------------------------------------------------------------
#isMonteCarlo = False		#if running on Monte-Carlo file
doDAvertex   = False		#if you want to enable the AnnealingVertexing + stores it
doMBTracking = False		#does the low-pt tracking + merging with generalTracks + vertices + stores them. If false, store only generalTracks & offlinePrimaryVertices.
useMITFilter = False		#Not the official CMS one, but cleans better. Needed for low-pt tracking. If off, no filter is used. //the standard "noscraping" filter is used.
storeJets    = True		#stores jets. Need to choose which collections and corrections to store, and if do the Dijets.
storeMET     = True		#stores "met", "pfMet" and "tcMet". If isMonteCarlo=true , also stores "genMetTrue".
storeCastor  = True 		#stores castorrechits , castorjets
keepCMSData  = False 		#make another CMSSW root files with all the collections from the input file and those created in the path

###
storePFCands     = True
storeCaloObjects = True
storeMuons       = True
storeElectrons   = True
storeZDC         = True
storeFSC         = False
