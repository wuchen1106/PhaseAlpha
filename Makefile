CXX = g++
BINARIESDIR = BinaryFiles
BINDIR = $(BINARIESDIR)/bin
APPDIR = app
SRCDIR = src
SUFFIX = .cxx
HEADER = .hxx
EXECUTE = 

LDFLAGS  += -shared
#CXXFLAGS += -W -Wall -O3 -DNDEBUG -fPIC
CXXFLAGS += -W -Wall -O2 -g -DNDEBUG -fPIC
CXXFLAGS += $(shell root-config --cflags)
LIBS     += $(shell $(ROOTSYS)/bin/root-config --glibs) -lMinuit
LIBS     += -pthread -lm -ldl -rdynamic -lGeom -lEG
INCS     += -I$(shell $(ROOTSYS)/bin/root-config --incdir)

SRCA := $(wildcard $(APPDIR)/*$(SUFFIX))
SRCS := $(wildcard $(SRCDIR)/*$(SUFFIX))
HEADERS := $(wildcard $(SRCDIR)/*$(HEADER))

TGTS = $(addprefix $(BINDIR)/, $(notdir $(basename $(SRCA))))

.PHONY: all
all: $(TGTS)

$(BINDIR)/%: $(APPDIR)/%$(SUFFIX) $(HEADERS)
	mkdir -p $(BINDIR); \
	$(CXX) $(CXXFLAGS) $(LIBS) $(INCS) $(APPDIR)/$(notdir $@)$(SUFFIX) -o $@${EXECUTE} $(filter-out $(APPDIR)/%$(SUFFIX), $^)

.PHONY: clean
clean:	
	rm -rf $(BINARIESDIR)
