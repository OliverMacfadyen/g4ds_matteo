#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.16...3.21)
project(g4ds11)
#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Find ROOT (required package)
#
#find_package(ROOT REQUIRED)
#IF (${ROOT_FOUND})
#   MESSAGE(STATUS "Found ROOT")
#ELSE (${ROOT_FOUND})
#   MESSAGE(STATUS "Could not locate ROOT")
#ENDIF (${ROOT_FOUND})

#----------------------------------------------------------------------------
# DEVELOPMENT MODE: WARNING = ERROR
#
set(ENV{G4DEV} 1)

IF ($ENV{G4DEV})
    MESSAGE(STATUS "Development mode: warning = error")
    add_compile_options(-Wall -Wextra -pedantic -Werror)
ENDIF()


#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include
                    ${ROOT_INCLUDE_DIRS})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)


#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(g4ds g4ds.cc ${sources} ${headers})
target_link_libraries(g4ds ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build B1. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(EXAMPLEG4DS_SCRIPTS
  macros/S1S2.mac
  macros/ambe.mac
  macros/amc.mac
  macros/argo.mac
  macros/bias.mac
  macros/calis.mac
  macros/cosmicray.mac
  macros/dart-external.mac
  macros/dev.mac
  macros/ds10.mac
  macros/ds20fluka.mac
  macros/ds20k_veto_optics.mac
  macros/ds20kbg.mac
  macros/ds20klm.mac
  macros/ds20klsv.mac
  macros/ds50lsv.mac
  macros/dsProto.mac
  macros/dsProtoProto.mac
  macros/gantry.mac
  macros/geo.mac
  macros/hepevt.mac
  macros/licorne.mac
  macros/multi.mac
  macros/multisp.mac
  macros/neutrons.mac
  macros/neutron_captures.mac
  macros/rdm.mac
  macros/rdmchain.mac
  macros/red.mac
  macros/red_lns.mac
  macros/run.mac
  macros/scs.mac
  macros/vis.mac
  macros/viswrl.mac
  macros/writeGDML.mac
  )

foreach(_script ${EXAMPLEG4DS_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS g4ds DESTINATION bin)
