# SVN: $Id: Makefile 243 2013-08-27 01:44:59Z mushkar $

include $(ROOTSYS)/etc/Makefile.arch

IncSuf = h

FIT_SRC    = $(wildcard TFit*.$(SrcSuf)) FitDict.$(SrcSuf)
FIT_INC    = $(patsubst %.$(SrcSuf), %.$(IncSuf), $(wildcard TFit*.$(SrcSuf))) FitTypes.h FunComptonTheory.h
FIT_OBJS   = $(patsubst %.$(SrcSuf), ../obj/%.$(ObjSuf), $(FIT_SRC))
FIT_SO     = ../lib/libFit.$(DllSuf)
FIT_LIB    = $(shell pwd)/$(FIT_SO)
INCFLAGS     = -I$(ROOTSYS)/include

# PROGS
FIT_MAINO  = ../obj/FitMain.$(ObjSuf)
FIT_BIN    = ../bin/Fit
#
PartDeriv_MAINO  = ../obj/PartDerivMain.$(ObjSuf)
PartDeriv_BIN    = ../bin/PartDeriv

LIBS += -lMathMore
# LDFLAGS += -fopenmp
CXXFLAGS += -g
# CXXFLAGS += -g -fopenmp

all: Fit PartDeriv docs

Fit: $(FIT_BIN)
	@echo "$@ done"
	@echo "##-----------------------------------------------------------"

PartDeriv: $(PartDeriv_BIN)
	@echo "$@ done"
	@echo "##-----------------------------------------------------------"

clean:
	@rm -fv $(FIT_OBJS)
	@rm -fv $(FIT_MAINO)
	@rm -fv $(FIT_BIN)
	@rm -fv $(PartDeriv_MAINO)
	@rm -fv $(PartDeriv_BIN)
	@rm -fv $(FIT_SO)
	@rm -fv FitDict.$(SrcSuf) FitDict.$(IncSuf)
	@rm -rfv $(Fit)/htmldoc

docs:
	@echo "Invoking ROOT html automatic documentation..."
# 	@rm -rf $(Fit)/htmldoc
	root -b -n -q -l MakeHtmlDoc.C
	@echo "$@ done"
	@echo "##-----------------------------------------------------------"

print:
	@echo $(FIT_SRC)
	@echo "##-----------------------------------------------------------"
	@echo $(wildcard TFit*.$(SrcSuf))
	@echo "##-----------------------------------------------------------"
	@echo $(FIT_INC)
	@echo "##-----------------------------------------------------------"
	@echo $(FIT_OBJS)
	@echo "$@ done"
	@echo "##-----------------------------------------------------------"

$(FIT_SO): $(FIT_OBJS)
	$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(OutPutOpt) $@ $(EXPLLINKLIBS)
	@echo "$@ done"
	@echo "##-----------------------------------------------------------"


$(FIT_BIN): $(FIT_MAINO) $(FIT_SO)
	$(LD) $(LDFLAGS) $(FIT_MAINO) $(FIT_LIB) $(LIBS) $(OutPutOpt)$@
	$(MT_EXE)
	@echo "$@ done"
	@echo "##-----------------------------------------------------------"
	
$(PartDeriv_BIN): $(PartDeriv_MAINO) $(FIT_SO)
	$(LD) $(LDFLAGS) $(PartDeriv_MAINO) $(FIT_LIB) $(LIBS) $(OutPutOpt)$@
	$(MT_EXE)
	@echo "$@ done"
	@echo "##-----------------------------------------------------------"

FitDict.$(SrcSuf): $(FIT_INC) FitLinkDef.$(IncSuf)
	@echo "Generating dictionary $@..."
	@rootcint  -f $@ -c $(INCFLAGS) $^
	@echo "$@ done"
	@echo "##-----------------------------------------------------------"

../obj/%.$(ObjSuf): %.$(SrcSuf)
	$(CXX) $(CXXFLAGS) -o $@ -c $<
