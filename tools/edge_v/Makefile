##
# @file This file is part of EDGE.
#
# @author Junyi Qiu (juq005 AT ucsd.edu)
# @author Alexander Breuer (anbreuer AT ucsd.edu)
#
# @section LICENSE
# Copyright (c) 2017, Regents of the University of California
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# @section DESCRIPTION
# This is the Makefile of Edge-V.
##

CXX ?= mpicxx

PREFIX ?= .


-include ./Makefile.inc

UCVMINC=$(UCVMDIR)/include
UCVMLIB=$(UCVMDIR)/lib

PROJ4INC=$(PROJ4DIR)/include
PROJ4LIB=$(PROJ4DIR)/lib

UCVMCLDFLAGS=-L$(UCVMDIR)/model/cvms426/lib  -lcvmsi       -Wl,-rpath,$(UCVMDIR)/model/cvms426/lib  \
             -L$(UCVMDIR)/model/cvms/lib     -lcvms        -Wl,-rpath,$(UCVMDIR)/model/cvms/lib     \
             -L$(UCVMDIR)/model/cvmh1510/lib -lvxapi -lgeo -Wl,-rpath,$(UCVMDIR)/model/cvmh1510/lib \
             -L$(UCVMDIR)/model/cencal/lib   -lcencalvm    -Wl,-rpath,$(UCVMDIR)/model/cencal/lib   \
             -L$(UCVMDIR)/lib/euclid3/lib    -letree       -Wl,-rpath,$(UCVMDIR)/lib/euclid3/lib    \
             -L$(UCVMDIR)/lib/fftw/lib       -lfftw3       -Wl,-rpath,$(UCVMDIR)/lib/fftw/lib

MOABINC=$(MOABDIR)/include
MOABLIB=$(MOABDIR)/lib

-include $(MOABLIB)/moab.make
MOAB_LINKS := $(shell echo $(MOAB_LIBS_LINK) | sed 's/-all-static//g')

CXXFLAGS=-I. -I$(PROJ4INC) -I$(UCVMINC) -I$(MOABINC) -fopenmp -std=c++11 -Wall
LDFLAGS=-L$(PROJ4LIB) -lproj -Wl,-rpath,$(PROJ4LIB)                \
        -L$(UCVMLIB)  -lucvm $(UCVMCLDFLAGS) -Wl,-rpath,$(UCVMLIB) \
        $(MOAB_LINKS) -lm -ldl

SRCDIR=./src
OBJDIR=./build
RELEASE=$(PREFIX)/bin

_DEPS=vm_utility.h
_OBJ=vm_utility.o edge_v.o
_POS=vm_utility.o write_pos.o

DEPS=$(patsubst %,$(SRCDIR)/%,$(_DEPS))
OBJ=$(patsubst %,$(OBJDIR)/%,$(_OBJ))
POS=$(patsubst %,$(OBJDIR)/%,$(_POS))

all: build write_pos

build: $(OBJDIR)/edge_v

write_pos: $(OBJDIR)/write_pos

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(DEPS)
	if [ ! -d $(OBJDIR)  ]; then mkdir -p $(OBJDIR) ; fi
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(OBJDIR)/edge_v: $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
	if [ ! -d $(RELEASE) ]; then mkdir -p $(RELEASE); fi
	cp $(OBJDIR)/edge_v $(RELEASE)/

$(OBJDIR)/write_pos: $(POS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
	if [ ! -d $(RELEASE) ]; then mkdir -p $(RELEASE); fi
	cp $(OBJDIR)/write_pos $(RELEASE)/

.PHONY: cleanobj clean

cleanobj:
	rm -f $(OBJDIR)/*.o $(OBJDIR)/edge_v $(OBJDIR)/write_pos

clean:
	rm -f $(OBJDIR)/*.o $(OBJDIR)/edge_v $(OBJDIR)/write_pos $(RELEASE)/edge_v $(RELEASE)/write_pos
