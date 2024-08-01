# $Id: GNUmakefile,v 1.5 2015/10/01 08:46:12 dfranco Exp $
# --------------------------------------------------------------
# If you want to use the graphics (visualisation) in the programm
# just comment the line : export G4VIS_NONE=1
# and recompile the program by 2 commands:
# make & make clean
# ---------------------------------------------------------------

G4BIN := ./

name := g4ds

G4TARGET := $(name)
G4EXLIB := true

G4VIS_USE := 0
BX_USE_OPERA := 0

#ifndef G4INSTALL
  #G4INSTALL = ../../..
#G4INSTALL = $(GEANT4_FQ_DIR)/share/Geant4-10.0.0/geant4make/
#G4INSTALL = $(GEANT4_FQ_DIR)/share/Geant4-10.0.1/geant4make/
#endif

.PHONY: all
all: lib bin

HOST = $(shell hostname)
CNAF := $(shell hostname | grep 'cnaf')
LYON := $(shell hostname | grep 'cca')
ARGUS = argus.Princeton.EDU
DEATHSTAR = deathstar
YODA = phy-yoda.Princeton.EDU
CPPFLAGS += -I$(ROOTSYS)/include -DBX_USE_OPERA 
ifneq ($(LD_LIBRARY_PATH),) # linux
  EXTRALIBS = $(shell echo "-L$(LD_LIBRARY_PATH)" | sed 's/:/ -L/g')
endif

ifeq ($(HOST),$(ARGUS))
EXTRALIBS = $(shell root-config --glibs) /sw/lib/libxerces-c.27.dylib
else
ifeq ($(HOST),$(DEATHSTAR))
EXTRALIBS = -lxerces-c -L/usr/local/bin/root5.34/lib -lMathCore -lCore -lCint -lHist -lNet -lGraf -lGraf3d -lGpad -lTree -lRint -lMatrix -lPhysics -lThread -lGui -pthread -rdynamic
CPPFLAGS += -I/usr/local/root_v5.32.04/include/root/
else
ifeq ($(HOST),$(YODA))
EXTRALIBS = $(shell root-config --glibs --cflags) -L/usr/local/root_v5.32.04/lib/root
CPPFLAGS += -I/usr/local/root_v5.32.04/include/root
else
ifeq ($(HOST),$(CNAF))
EXTRALIBS = $(shell root-config --glibs)
else
ifeq ($(HOST),$(LYON))
EXTRALIBS = $(shell root-config --glibs)
else
EXTRALIBS += $(shell root-config --glibs)
CPPFLAGS += $(shell root-config --cflags) -Wno-shadow -ggdb
endif
endif
endif
endif
endif



# type `make quiet=yes` to disable all warnings
ifeq ($(quiet),yes)
  CPPFLAGS+=-w
endif

include $(G4INSTALL)/config/binmake.gmk

visclean:
	rm -f g4*.prim g4*.eps
	rm -f ./Darwin-g++/*.wrl
	rm -f .DAWN_*

logclean:
	rm -f ./Linux-g++/*.fil
	rm -f ./Linux-g++/*.log
	rm -f ./Linux-g++/*.out
	rm -f ./Linux-g++/*.root
	rm -f ./Darwin-g++/*.fil
	rm -f ./Darwin-g++/*.log
	rm -f ./Darwin-g++/*.out
	rm -f ./Darwin-g++/*.root
