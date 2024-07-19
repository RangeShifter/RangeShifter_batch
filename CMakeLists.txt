cmake_minimum_required(VERSION 3.10)

# set the project name and version
project(RangeShifter VERSION 2.1.0)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# add RScore to source search path
set(batchmode ON) # don't config RScore as standalone
add_subdirectory(src/RScore)

# add the executable
add_executable(RangeShifter src/Main.cpp src/BatchMode.cpp)

target_compile_definitions(RangeShifter PRIVATE RSWIN64)

# enable LINUX_CLUSTER macro on Linux + macOS
if(CMAKE_SYSTEM_NAME STREQUAL "Linux" OR CMAKE_SYSTEM_NAME STREQUAL "Darwin")
	add_compile_definitions(LINUX_CLUSTER)
endif()

# Debug Mode by default, unless "release" is passed i.e. `cmake -Drelease=`
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

# link RScore to the executable as a library
target_link_libraries(RangeShifter PUBLIC RScore)

# add RScore as an include dir
target_include_directories(RangeShifter PUBLIC "${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}/RScore")
