IDIR       = include
ODIR       = build
SDIR       = Analyzer/bin

BASEDIR  = $(CMSSW_BASE)/src/SUSYBSMAnalysis/
BINDIR   = $(BASEDIR)/Analyzer/bin

CXX        = g++

#CXXFLAGS  += -I$(IDIR) -I$(CMSSW_BASE)/$(SDIR) -std=c++0x
CXXFLAGS  += -I$(IDIR) -I$(CMSSW_BASE)/$(SDIR) -std=c++0x
CXXFLAGS  += -g -O3 ## Optimization flag
CXXFLAGS  += $(shell root-config --cflags) ## Include ROOT

CXXDEPFLAGS = -MMD -MP

LD         = g++

LIBS       = $(shell root-config --glibs)

PROGRAMS = RunBackgroundPrediction


all: mkobj $(PROGRAMS)

mkobj:
	@mkdir -p $(ODIR)

$(ODIR)/%.o : $(SDIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -o $@ -c $<


RunBackgroundPrediction: $(ODIR)/BackgroundPrediction.o
	$(LD) $^ $(LIBS) -o $@

clean:
	rm -f $(ODIR)/*.o $(ODIR)/*.so $(ODIR)/*.d $(PROGRAMS)
