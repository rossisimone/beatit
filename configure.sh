# BUILD MODE
BUILD=RELEASE

# Source directory of BeatIt. Use full Path
BeatIt_SOURCE_DIR=

# Install directory
PREFIX=$BeatIt_SOURCE_DIR/../install

# MPI
MPICC=mpicc
MPICXX=mpicxx

# If you have one, this can be useful
# TPL Base Folder
# Example
# TPL=~/TPL
TPL=

# petsc
# Example: 
# PETSC=$TPL/petsc/3.7.4/opt/
PETSC=
PETSC_LIB=

# LibMesh
# Example: 
# LIBMESH=$TPL/libmesh/1.0.0/opt/
LIBMESH=

# VTK
# Example
# VTK=$TPL/vtk/7.0.0/lib/cmake/vtk-7.0
VTK=

# Boost
# Example
# BOOST=$TPL/boost/1.61.0/
BOOST=

# OPTIONS
VERBOSITY=ON

rm -rf CMake*
cmake \
  -D CMAKE_INSTALL_PREFIX=$PREFIX \
  -D CMAKE_BUILD_TYPE=$BUILD \
  -D BUILD_SHARED_LIBS=ON \
  -D CMAKE_C_COMPILER=$MPICC \
  -D CMAKE_CXX_COMPILER=$MPICXX \
  -D CMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} -std=c++11" \
  -D LIBMESH_DIR=$LIBMESH  \
  -D PETSC_LIBRARIES=$PETSC_LIB  \
  -D VTK_DIR=$VTK  \
  -D PETSC_DIR=$PETSC \
  -D BOOST_ROOT=$BOOST \
  -D CMAKE_VERBOSE_MAKEFILE=$VERBOSITY \
  $BeatIt_SOURCE_DIR






