TotemSrc=../TOTEMdataFormat/src
CMSSrc=../CMSdataFormat/src
TotemLibrary=../TOTEMdataFormat/lib/libTOTEMdataFormat.so
CMSLibrary=../CMSdataFormat/lib/libCMSdataFormat.so

CXXFLAGS+=$(shell root-config --cflags) -Wall -I. -I$(TotemSrc) -I$(CMSSrc)
LDFLAGS+=$(shell root-config --libs) -lTreePlayer 

all: ../TOTEMdataFormat/utilities/Makefile ../CMSdataFormat/utilities/Makefile
	cd ../TOTEMdataFormat/utilities/; make
	cd ../CMSdataFormat/utilities/; make
	make all2

all2: merge12 combine mergeNtuples findOrbitOffset mergeMCNtuples

merge12: utilities.h $(CMSLibrary) $(TotemLibrary)

combine: utilities.h $(CMSLibrary) $(TotemLibrary)

mergeNtuples: $(CMSLibrary) $(TotemLibrary)

mergeMCNtuples: $(CMSLibrary) $(TotemLibrary)

findOrbitOffset: $(CMSLibrary) $(TotemLibrary)

clean:
	cd ../TOTEMdataFormat/utilities/; make clean
	cd ../CMSdataFormat/utilities/; make clean
	rm -f merge12 combine mergeNtuples findOrbitOffset *.o

