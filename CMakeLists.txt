cmake_minimum_required(VERSION 3.16)

# set the project name and version
project(RangeShifter VERSION 2.1.0)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
option(BUILD_SHARED_LIBS "Build shared libraries" OFF)

# get SFML (graphics library) from GitHub
include(FetchContent)
FetchContent_Declare(SFML
    GIT_REPOSITORY https://github.com/SFML/SFML.git
    GIT_TAG 2.6.x)
FetchContent_MakeAvailable(SFML)

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
# link RScore to the executable as a library
target_link_libraries(RangeShifter PUBLIC RScore PRIVATE sfml-graphics)

# add RScore as an include dir
target_include_directories(RangeShifter PUBLIC "${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}/RScore")


# link SFML-graphics

#if(WIN32) # Windows
#    add_custom_command(
#        TARGET RangeShifter
#        COMMENT "Copy OpenAL DLL"
#        PRE_BUILD COMMAND ${CMAKE_COMMAND} -E copy ${SFML_SOURCE_DIR}/extlibs/bin/$<IF:$<EQUAL:${CMAKE_SIZEOF_VOID_P},8>,x64,x86>/openal32.dll $<TARGET_FILE_DIR:RangeShifter>
#        VERBATIM)
#endif()