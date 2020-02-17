cmake_minimum_required(VERSION 3.10)

### Hunter package manager ###
include("cmake/HunterGate.cmake")
HunterGate(
    URL "https://github.com/ruslo/hunter/archive/v0.23.153.tar.gz"
    SHA1 "d4c16d7474a1744eff33b832a82456852232ea60"
)

project(cc_codes)

# Check that this is not an in-source build:
if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
   message(SEND_ERROR "In-source builds are not allowed.")
endif()

#setup compiler
set(CMAKE_CXX_STANDARD 17)

#enable warnings
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG(-Wall COMPILER_SUPPORTS_WALL)
CHECK_CXX_COMPILER_FLAG(-Wextra COMPILER_SUPPORTS_WEXTRA)
CHECK_CXX_COMPILER_FLAG(-Wconversion COMPILER_SUPPORTS_WCONVERSION)

add_library(warning_flags INTERFACE)
ADD_COMPILE_OPTIONS(-Wpedantic -O3)
if(COMPILER_SUPPORTS_WALL)
  target_compile_options(warning_flags INTERFACE -Wall)
endif()
if(COMPILER_SUPPORTS_WEXTRA)
  target_compile_options(warning_flags INTERFACE -Wextra)
endif()
if(COMPILER_SUPPORTS_WCONVERSION)
  target_compile_options(warning_flags INTERFACE -Wconversion)
endif()

# Get Dependencies
# (don't forget to update cmake/Config.cmake.in !
###############################################################################
# Get Eigen
hunter_add_package(Eigen)
find_package(Eigen3 CONFIG REQUIRED)

# Get LehrFEMpp
hunter_add_package(lehrfempp)
find_package(lehrfempp CONFIG REQUIRED)

# Get Boost
hunter_add_package(Boost COMPONENTS chrono timer system filesystem program_options)
find_package(Boost CONFIG REQUIRED chrono timer system filesystem program_options)


# Include Subdirectories
###############################################################################
add_subdirectory(codes)
