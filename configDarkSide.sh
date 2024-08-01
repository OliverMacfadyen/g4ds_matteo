#!/bin/sh

#export G4VIS_USE_OPENGLX=1

#unset ROOTSYS
export G4DS=${PWD}

source /sps/nusol/software/miniconda3/envs/darkside/share/Geant4-11.0.2/geant4make/geant4make.sh


#export ROOTSYS=/usr/local/root/v5.34.11/

export DSDATA=$PWD/data/physics

export PATH=${G4DS}/tools:${PATH}
#export PATH=${ROOTSYS}/bin:${PWD}/tools:${PATH}
#export LD_LIBRARY_PATH=${ROOTSYS}/lib:${ROOTSYS}/lib/root:${LD_LIBRARY_PATH}

#source /usr/local/geant4/geant4.10.00.p01/share/Geant4-10.0.1/geant4make/geant4make.sh > /dev/null
#source /usr/local/geant4/geant4.9.6.p01/share/Geant4-9.6.1/geant4make/geant4make.sh

#export G4VRMLFILE_MAX_FILE_NUM=1
#export G4VIS_USE=1
#unset G4VIS_USE_DAWN           
#unset G4VIS_USE_OPENGLQT
#unset G4VIS_USE_RAYTRACERX
#unset G4VIS_USE_VRML
#unset G4VIS_USE_OIX
#unset G4VIS_USE_OPENGLXM
#unset G4UI_USE_QT

#export G4DATA=/usr/local/geant4/share

#export G4SAIDXSDATA=$G4DATA/G4SAIDDATA1.1
#export G4SAIDDATA=$G4DATA/G4SAIDDATA1.1
#export G4LEDATA=$G4DATA/G4EMLOW6.35
#export G4NEUTRONHPDATA=$G4DATA/G4NDL4.4
#export G4NEUTRONXSDATA=$G4DATA/G4NEUTRONXS1.4 
#export G4LEVELGAMMADATA=$G4DATA/PhotonEvaporation3.1
#export G4RADIOACTIVEDATA=$G4DATA/RadioactiveDecay4.2
#export NeutronHPCrossSections=$G4DATA/G4NDL4.4
#export G4REALSURFACEDATA=$G4DATA/RealSurface1.0
#export G4PIIDATA=$G4DATA/G4PII1.3

export G4WORKDIR=$PWD/.g4ds


#export LD_LIBRARY_PATH=$G4WORKDIR/tmp/Linux-g++/g4ds/:$XERCESCROOT/lib:$LD_LIBRARY_PATH
						

