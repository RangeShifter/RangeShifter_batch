# Config file for compilation with CMake

if (NOT batchmode) # that is, RScore as a standalone
	cmake_minimum_required(VERSION 3.10)
	# set the project name and version
	project(RScore VERSION 2.1.0)
	# specify the C++ standard
	set(CMAKE_CXX_STANDARD 17)
	set(CMAKE_CXX_STANDARD_REQUIRED True)
	add_executable(RScore Main.cpp Species.cpp Cell.cpp Community.cpp FractalGenerator.cpp Genome.cpp Individual.cpp Landscape.cpp Model.cpp Parameters.cpp Patch.cpp Population.cpp RandomCheck.cpp RSrandom.cpp SubCommunity.cpp Utils.cpp)
else() # that is, RScore compiled as library within RangeShifter_batch
	add_library(RScore Species.cpp Cell.cpp Community.cpp FractalGenerator.cpp Genome.cpp Individual.cpp Landscape.cpp Model.cpp Parameters.cpp Patch.cpp Population.cpp RandomCheck.cpp RSrandom.cpp SubCommunity.cpp Utils.cpp)
endif()

# pass config definitions to compiler
target_compile_definitions(RScore PRIVATE RSWIN64)

# enable LINUX_CLUSTER macro on Linux + macOS
if(CMAKE_SYSTEM_NAME STREQUAL "Linux" OR CMAKE_SYSTEM_NAME STREQUAL "Darwin")
	add_compile_definitions("LINUX_CLUSTER")
endif()

# Debug Mode by default, unless "release" is passed
if(NOT DEFINED release)
	add_compile_definitions(RSDEBUG)
endif()


if(DEFINED RS_CONTAIN)
	add_compile_definitions(RS_CONTAIN)
endif()

if(DEFINED TEMPMORT)
	add_compile_definitions(TEMPMORT)
endif()

if(DEFINED GROUPDISP)
	add_compile_definitions(GROUPDISP)
endif()

if(DEFINED PARTMIGRN)
	add_compile_definitions(PARTMIGRN)
endif()

if(DEFINED SEASONAL)
	add_compile_definitions(SEASONAL)
endif()

if(DEFINED RS_ABC)
	add_compile_definitions(RS_ABC)
endif()

if(DEFINED GOBYMODEL)
	add_compile_definitions(GOBYMODEL)
endif()

if(DEFINED VIRTUALECOLOGIST)
	add_compile_definitions(VIRTUALECOLOGIST)
endif()

if(DEFINED BUTTERFLYDISP)
	add_compile_definitions(BUTTERFLYDISP)
endif()

if(DEFINED PEDIGREE)
	add_compile_definitions(PEDIGREE)
endif()

if(NOT batchmode)
	target_include_directories(RScore PUBLIC "${PROJECT_BINARY_DIR}")
endif() 