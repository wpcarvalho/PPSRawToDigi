ROOT_INCS = -I$(ROOTSYS)/include
INCS = -I $(CMSSW_RELEASE_BASE)/src -I$(CMSSW_BASE)/src $(ROOT_INCS)
INCS += -I/afs/cern.ch/cms/sw/slc4_ia32_gcc345/external/clhep/1.9.3.2/include/
INCS +=	-I/afs/cern.ch/cms/sw/slc4_ia32_gcc345/external/boost/1.34.1/include/
INCS += -I/afs/cern.ch/cms/sw/slc4_ia32_gcc345/external/hepmc/2.01.10/include/

#ROOT_LIBS  = -L$(ROOTSYS)/lib -lCore -lCint -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lz -pthread -lm -ldl -rdynamic -lGui
ROOT_LIBS  = `$(ROOTSYS)/bin/root-config --libs --glibs`
LIBS  = $(ROOT_LIBS)
LIBS +=	-L$(CMSSW_BASE)/lib/$(SCRAM_ARCH) -L$(CMSSW_RELEASE_BASE)/lib/$(SCRAM_ARCH)
LIBS +=	-L/afs/cern.ch/cms/sw/slc4_ia32_gcc345/external/zlib/1.1.4/lib/
LIBS +=	-lSimDataFormatsHepMCProduct -lSimDataFormatsHepMCProductCapabilities
LIBS +=	-lFWCoreFWLite
LIBS +=	-lRecoTotemRPRPRecoDataFormats
LIBS +=	-lDataFormatsTotemRPDetId

# the following line runs colorizer wrapper around /usr/cmssw/slc4_ia32_gcc345/external/gcc-wrapper/3.4.5m32/bin/c++
#CXX = gcc_cmssw
#CXX = gcc
CXX = /usr/cmssw/slc4_ia32_gcc345/external/gcc-wrapper/3.4.5m32/bin/c++

all: runAnalysis prob

runAnalysis: analyze.o
	@echo LINKING analyze
	@$(CXX) $(LIBS) analyze.o -orunAnalysis

prob: prob.o
	@echo LINKING prob
	@$(CXX) $(LIBS) prob.o -oprob

hough: hough.cc
	@$(CXX) $(ROOT_INCS) $(ROOT_LIBS) hough.cc -ohough

%.o: %.cc
	@echo COMPILING $@
	@$(CXX) $(INCS) -c $< -o $@
