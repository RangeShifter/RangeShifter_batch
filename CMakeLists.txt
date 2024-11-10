cmake_minimum_required(VERSION 3.16)

# set the project name and version
project(RangeShifter VERSION 2.1.0)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

option(BUILD_SHARED_LIBS "Build shared libraries" OFF) # required by SFML's own build

# get SFML (graphics library) from GitHub
include(FetchContent)
FetchContent_Declare(SFML
    GIT_REPOSITORY https://github.com/SFML/SFML.git
    GIT_TAG 2.6.x)
FetchContent_MakeAvailable(SFML)

# add RScore to source search path
set(batchmode ON) # don't config RScore as standalone
set(batchview ON) # tell RScore to enable BATCH_VIEW macro
add_subdirectory(src/RScore)

# declare target executable
add_executable(RangeShifter src/Main.cpp src/BatchMode.cpp src/BatchView.cpp)

# enable LINUX_CLUSTER macro on Linux + macOS
if(CMAKE_SYSTEM_NAME STREQUAL "Linux" OR CMAKE_SYSTEM_NAME STREQUAL "Darwin")
	add_compile_definitions(LINUX_CLUSTER)
endif()

# Debug Mode by default, unless "release" is passed i.e. `cmake -Drelease=`
if(NOT DEFINED release)
	add_compile_definitions(RSDEBUG)
endif()

# link external dependencies
target_link_libraries(RangeShifter 
						PUBLIC RScore # core code for the simulation
						PRIVATE sfml-graphics) # graphics library

# add RScore as an include dir
target_include_directories(RangeShifter PUBLIC "${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}/RScore")