CMAKE_MINIMUM_REQUIRED(VERSION 3.10)

project(one_d) 

set(CMAKE_CXX_STANDARD 17)

# find the eigen package
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

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

# define executables
add_executable(conv src/conv_main.cc
		    src/galerkin_mat.cc)
add_executable(main src/main.cc
	            src/galerkin_mat.cc)

target_link_libraries (main Eigen3::Eigen)
target_link_libraries (conv Eigen3::Eigen)
