cmake_minimum_required(VERSION 3.10)

# set the project name and version
project(RangeShifter VERSION 2.1.0)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# add RScore to source search path
add_subdirectory(src/RScore)

# add the executable
add_executable(RangeShifter src/main.cpp src/BatchMode.cpp)

# pass config definitions to compiler
target_compile_definitions(RangeShifter PRIVATE "RSDEBUG" "RSWIN64" "BATCH")

# link RScore to the executable as a library
target_link_libraries(RangeShifter PUBLIC RScore)

# add RScore as an include dir
target_include_directories(RangeShifter PUBLIC "${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}/RScore")