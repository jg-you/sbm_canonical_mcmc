# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Project
# ~~~~~~~~~~~~~~~~~~~~~~~~~
cmake_minimum_required(VERSION 2.8)
project(sbm_inference)

# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Check for functionalities 
# ~~~~~~~~~~~~~~~~~~~~~~~~~

# Check for C++ 11 with CMAKE v.2.8 
# Note: there exists a cleaner method in CMAKE 3.1+ [target_compile_features()]
#       but stable Debian distributions tend to distribute CMAKE 2.8 as the default.
#
# Solution found @ http://stackoverflow.com/a/25836953/1851837.
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
    message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

# Boost
find_package( Boost 1.40 COMPONENTS program_options )
if (Boost_FOUND)
    set(HAVE_LIBBOOST_PROGRAM_OPTIONS 1)
else()
    message(STATUS "C++ Boost::program_options is not installed. Using limited interface.")
    set(HAVE_LIBBOOST_PROGRAM_OPTIONS 0)
endif ()

# Steady clock (Google code)
include(cmake_tests/CXXFeatureCheck.cmake)
# If successful, then HAVE_STEADY_CLOCK is set to 1
set(HAVE_STEADY_CLOCK 0)
cxx_feature_check(STEADY_CLOCK)

# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Options
# ~~~~~~~~~~~~~~~~~~~~~~~~~
option(LOGGING 
        "Log input information to std::clog." ON)

# Defaults
set (LOGGING 1)
if (!LOGGING)
  set (LOGGING 0)
endif()

# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Build
# ~~~~~~~~~~~~~~~~~~~~~~~~~

configure_file (
  "${PROJECT_SOURCE_DIR}/src/config.h.in"
  "${PROJECT_BINARY_DIR}/src/config.h"
  )
include_directories("${PROJECT_BINARY_DIR}")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_subdirectory(src)

