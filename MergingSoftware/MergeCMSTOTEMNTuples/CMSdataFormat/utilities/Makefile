SRC=../src
OBJ=../lib

CXXFLAGS+=$(shell root-config --cflags) -Wall -I$(SRC) -ggdb3 -fPIC
LDFLAGS+=$(shell root-config --glibs)  -ggdb3 

HEAD=$(wildcard $(SRC)/My*.h)
OBJS1=$(HEAD:.h=.o)
OBJS=$(subst $(SRC),$(OBJ),$(OBJS1))

LINKD=$(SRC)/LinkDef


all: $(OBJ) $(OBJ)/libCMSdataFormat.so

$(OBJ)/libCMSdataFormat.so: $(OBJS) $(OBJ)/eventdict.o 
	$(CXX) $(LDFLAGS) -shared -o $@ $^	
	@echo ""
	#to load the library to CINT one needs X connection
	@echo ""

$(OBJ)/%.o: $(SRC)/%.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@
	@echo ""

$(SRC)/eventdict.cc: $(HEAD) $(LINKD).h
	rootcint -f $@ -c -p $^
	@echo ""
	#compatibility with cmssw
	mv $(LINKD).h $(LINKD)h
	@echo ""

$(LINKD).h:
	#compatibility with cmssw
	mv $(LINKD)h  $(LINKD).h
	@echo ""
$(OBJ):
	mkdir -p $(OBJ)
clean:
	rm -rf $(OBJ)

