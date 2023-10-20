cmake_minimum_required(VERSION 3.10)

# set the project name and version
project(RangeShifter VERSION 2.1)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# configure a header file to pass some of the CMake settings
# to the source code
#configure_file(TutorialConfig.h.in TutorialConfig.h)

# add RScore to source search path
add_subdirectory(src/RScore)

# add the executable
add_executable(RangeShifter src/main.cpp src/BatchMode.cpp src/Version.cpp)

# link RScore to the executable as a library
target_link_libraries(RangeShifter PUBLIC RScore)

# add RScore as an include dir
target_include_directories(RangeShifter PUBLIC "${PROJECT_SOURCE_DIR}/RScore")