cmake_minimum_required(VERSION 3.10)

# set the project name and version
project(RangeShifter VERSION 2.1.0)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# add RScore to source search path
set(batchmode ON) # don't config RScore as standalone

if(DEFINED OMP)
	set(OMP ON)
endif()

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

# Enable parallelisation with OpenMP
# define macro during config: `cmake -DOMP=`
if(DEFINED OMP)
	find_package(OpenMP COMPONENTS CXX)
	if(OpenMP_CXX_FOUND)
		target_link_libraries(RangeShifter PUBLIC OpenMP::OpenMP_CXX)
	endif()
endif()

# link RScore to the executable as a library
target_link_libraries(RangeShifter PUBLIC RScore)

# add RScore as an include dir
target_include_directories(RangeShifter PUBLIC "${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}/RScore")
