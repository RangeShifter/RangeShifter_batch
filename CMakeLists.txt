cmake_minimum_required(VERSION 3.10)

# set the project name and version
project(RangeShifter VERSION 2.1)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# configure a header file to pass some of the CMake settings
# to the source code
#configure_file(TutorialConfig.h.in TutorialConfig.h)

# add_subdirectory(src/RScore)

# add the executable
add_executable(RangeShifter src/main.cpp src/BatchMode.cpp src/Version.cpp src/RScore/Species.cpp src/RScore/Cell.cpp src/RScore/Community.cpp src/RScore/FractalGenerator.cpp src/RScore/ src/RScore/Genome.cpp src/RScore/Individual.cpp src/RScore/Landscape.cpp src/RScore/Model.cpp src/RScore/Parameters.cpp src/RScore/Patch.cpp src/RScore/Population.cpp src/RScore/RandomCheck.cpp src/RScore/RSrandom.cpp src/RScore/SubCommunity.cpp)

# TODO 3: Use target_link_libraries to link the library to our executable
#target_link_libraries(RangeShifter PUBLIC MathFunctions)

# TODO 4: Add MathFunctions to Tutorial's target_include_directories()
# Hint: ${PROJECT_SOURCE_DIR} is a path to the project source. AKA This folder!
# add the binary tree to the search path for include files
# so that we will find TutorialConfig.h
# target_include_directories(RangeShifter PUBLIC "${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}/MathFunctions")