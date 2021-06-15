
#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(TE)

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
# ROOT SHIT ( see https://root-forum.cern.ch/t/using-cmake-to-build-geant4-root-project/20624/2 )
# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS RIO Net)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

include_directories(${ROOT_INCLUDE_DIRS})


#----------------------------------------------------------------------------
# Add option to debug
#
option(DEBUG_MAGNETIC_FIELD "Debug magnetic setting" OFF)
if(DEBUG_MAGNETIC_FIELD)
  add_definitions(-DDEBUG_MAGNETIC_FIELD)
endif()

option(DEBUG_ELECTRIC_FIELD "Debug electric setting" OFF)
if(DEBUG_ELECTRIC_FIELD)
  add_definitions(-DDEBUG_ELECTRIC_FIELD)
endif()


#----------------------------------------------------------------------------
# Add option to select particles from ROOT file for GEANT4 use
# To create ROOT file with fixed starting positions/energies/momenta
# compile and run generate_electrons.cc
#
option(GEANT4_ROOT "GEANT4 use" OFF)
if(GEANT4_ROOT)
  add_definitions(-DGEANT4_ROOT)
endif()

#----------------------------------------------------------------------------
# Add option to fire electrons for the iso_curves generation
# Energies are uniformly sampled between a minimum and a maximum (default: 1-150eV, dE = 1eV)
#
option(ISO_CURVES_GENERATOR "iso_curves generation" ON)
if(ISO_CURVES_GENERATOR)
  add_definitions(-DISO_CURVES_GENERATOR)
endif()

#----------------------------------------------------------------------------
# Add option to fire electrons for the iso_curves generation
# A finite range of energies is sampled around our standard energies (default ΔE = 5, centered on each standardE)
#
option(ENERGIES_RANGE_ISO_CURVES_GENERATOR "iso_curves generation for energy ranges around standard KE's" OFF)
if(ENERGIES_RANGE_ISO_CURVES_GENERATOR)
  add_definitions(-DENERGIES_RANGE_ISO_CURVES_GENERATOR)
endif()


#----------------------------------------------------------------------------
# Add option to select particles from ROOT file
# To create ROOT file with fixed starting positions/energies/momenta
# compile and run generate_electrons.cc
#
option(FIXED_PARTICLES "Particles from ROOT file" OFF)
if(FIXED_PARTICLES)
  add_definitions(-DFIXED_PARTICLES)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(test_electron test_electron.cc ${sources} ${headers})
target_link_libraries(test_electron ${ROOT_LIBRARIES})
target_link_libraries(test_electron ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build TE. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(EXAMPLETE_SCRIPTS
#   exampleTE.in
  gui.mac
  init_vis.mac
  vis.mac
  repeat.mac
  run_once.mac
  10k_repeat.mac
#  base.txt
  )

foreach(_script ${EXAMPLETE_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS test_electron DESTINATION bin)

