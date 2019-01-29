# minimum version of CMake is 3.3
cmake_minimum_required(VERSION 3.3 FATAL_ERROR)

# define project
project( gamma 
	 VERSION 0.9
	 LANGUAGES C CXX)

# we need at least C++11
set (CMAKE_CXX_STANDARD 11)

# CTF gives annoying format warnings
set (CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wno-format")

# Tell it where to find libs
set (LIBPATH "/usr/local")
set (INCPATH "/usr/local")

# required for unit testing
include(CTest)
enable_testing()

# required modules
include(GNUInstallDirs)
include(CheckIncludeFile)
include(CheckIncludeFileCXX)
include(CheckIncludeFiles)
include(ExternalProject)

find_package(Eigen3 3.3 REQUIRED)
find_package(Threads REQUIRED)
find_package(MPI REQUIRED)

    find_package(BLAS REQUIRED)
    include_directories(${BLAS_INCLUDE_PATH})


include_directories(${MPI_INCLUDE_PATH})

if (OPENMP_FLAGS)
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OPENMP_FLAGS}")
endif()

# build external projects
include(external/external.cmake)

# include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include/gamma)

# add subdirectories
add_subdirectory(src)
add_subdirectory(test)

# installation
install(DIRECTORY share/ DESTINATION share)
