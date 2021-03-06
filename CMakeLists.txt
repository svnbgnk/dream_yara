# ===========================================================================
#          DREAM-Yara - a distributed version of the Yara read mapper
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)
string(ASCII 27 Esc)
set(ColourBold  "${Esc}[1m")
set(ColourReset "${Esc}[m")
set(ColourRed   "${Esc}[31m")

message (STATUS "${ColourBold}Compiler Detection...${ColourReset}")

project (DREAM-YARA CXX)
# ----------------------------------------------------------------------------
# Make "Release" the default cmake build type
# ----------------------------------------------------------------------------

if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release CACHE STRING
    "Choose the type of build, options are: Debug Release RelWithDebInfo"
    FORCE)
endif ()

# ----------------------------------------------------------------------------
# Begin of dependency detection
# ----------------------------------------------------------------------------
message ("\n${ColourBold}Dependency detection${ColourReset}")

if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/include/seqan/include/seqan/version.h")
    set (CMAKE_INCLUDE_PATH
        ${CMAKE_CURRENT_SOURCE_DIR}/include/seqan/include
        ${CMAKE_INCLUDE_PATH})
    set (CMAKE_PREFIX_PATH
        ${CMAKE_CURRENT_SOURCE_DIR}/include/seqan/util/cmake
        ${CMAKE_PREFIX_PATH})
    set (CMAKE_MODULE_PATH
       ${CMAKE_CURRENT_SOURCE_DIR}/include/seqan/util/cmake
       ${CMAKE_MODULE_PATH})
    message (STATUS "Found a local SeqAn library provided with the DREAM-YARA source code.")
    message (    "   This will be preferred over system global headers.")
endif ()

if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/include/sdsl-lite/include/sdsl/bits.hpp")
    set (SDSL_INCLUDE_DIRS
        ${CMAKE_CURRENT_SOURCE_DIR}/include/sdsl-lite/include)
    message (STATUS "Found a local SDSL library provided with the DREAM-YARA source code.")
endif ()



# ----------------------------------------------------------------------------
# Add Yara targets
# ----------------------------------------------------------------------------

message (STATUS "${ColourBold}Configuring DREAM-Yara...${ColourReset}")
add_subdirectory(src)

# ----------------------------------------------------------------------------
# Warn if cmake build type is not "Release"
# ----------------------------------------------------------------------------

if (NOT CMAKE_BUILD_TYPE STREQUAL Release)
  message (STATUS "${ColourRed}CMAKE_BUILD_TYPE is not \"Release\", your binaries will be slow.${ColourReset}")
endif ()

# ----------------------------------------------------------------------------
# Add Tests
# ----------------------------------------------------------------------------

# message (STATUS "${ColourBold}Configuring DREAM-YARA Tests...${ColourReset}")

enable_testing ()

if (MODEL MATCHES ".*MemCheck.*")
    set (_VALGRIND_FLAG --valgrind)
else ()
    set (_VALGRIND_FLAG)
endif ()
find_package (PythonInterp)
if (PYTHONINTERP_FOUND)
  add_test (NAME app_test_dream_yara
            COMMAND python2
                    ${CMAKE_SOURCE_DIR}/tests/run_tests.py
                    ${_VALGRIND_FLAG}
                    ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR})
endif (PYTHONINTERP_FOUND)