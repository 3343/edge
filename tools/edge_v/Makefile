CXX=mpicxx

PREFIX ?= .

include ./Makefile.inc

UCVMINC=$(UCVMDIR)/include
UCVMLIB=$(UCVMDIR)/lib

PROJ4INC=$(PROJ4DIR)/include
PROJ4LIB=$(PROJ4DIR)/lib

UCVMCLDFLAGS=-L$(UCVMDIR)/model/cvms426/lib -lcvmsi \
             -L$(UCVMDIR)/model/cvms/lib -lcvms \
             -L$(UCVMDIR)/model/cvmh1510/lib -lvxapi -lgeo \
             -L$(UCVMDIR)/model/cencal/lib -lcencalvm \
             -L$(UCVMDIR)/lib/euclid3/lib -letree \
             -L$(UCVMDIR)/lib/fftw/lib -lfftw3

MOABINC=$(MOABDIR)/include
MOABLIB=$(MOABDIR)/lib

include $(MOABLIB)/moab.make
MOAB_LINKS := $(shell echo $(MOAB_LIBS_LINK) | sed 's/-all-static//g')

CXXFLAGS=-I. -I$(PROJ4INC) -I$(UCVMINC) -I$(MOABINC) -fopenmp 
LDFLAGS=-L$(PROJ4LIB) -lproj -L$(UCVMLIB) -lucvm $(UCVMCLDFLAGS) $(MOAB_LINKS) -lm -ldl


SRCDIR=./src
OBJDIR=./build
RELEASE=$(PREFIX)/bin

_DEPS=vm_utility.h
_OBJ =vm_utility.o edge_v.o 

DEPS=$(patsubst %,$(SRCDIR)/%,$(_DEPS))
OBJ =$(patsubst %,$(OBJDIR)/%,$(_OBJ))

all: build

build: $(OBJDIR)/edge_v

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(DEPS)
	if [ ! -d $(OBJDIR)  ]; then mkdir -p $(OBJDIR) ; fi
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(OBJDIR)/edge_v: $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
	if [ ! -d $(RELEASE) ]; then mkdir -p $(RELEASE); fi
	cp $(OBJDIR)/edge_v $(RELEASE)/

.PHONY: directory cleanobj clean

directory:
	if [ ! -d $(OBJDIR)  ]; then mkdir -p $(OBJDIR) ; fi
	if [ ! -d $(RELEASE) ]; then mkdir -p $(RELEASE); fi


cleanobj:
	rm -f $(OBJDIR)/*.o $(OBJDIR)/edge_v

clean:
	rm -f $(OBJDIR)/*.o $(OBJDIR)/edge_v $(RELEASE)/edge_v
