# CMakeLists.txt — build the ESMF two-component coupling example
cmake_minimum_required(VERSION 3.18)
project(ESMF_TwoComponentExample LANGUAGES Fortran)

# --- User configuration section ------------------------------------------
# Set this to your ESMF installation path if not found automatically.
# For example:
#   cmake -DESMF_DIR=/path/to/esmf .
set(ESMF_DIR "/software" CACHE PATH "Path to ESMF installation directory")

# NetCDF
set(NETCDF_LIBRARY "/software/lib/libnetcdf.so")
set(NETCDFF_LIBRARY "/software/lib/libnetcdff.so")
# -------------------------------------------------------------------------
# Find MPI (required by ESMF)
find_package(MPI REQUIRED Fortran)

# Locate ESMF library and module paths
if (NOT ESMF_DIR)
  message(FATAL_ERROR "Please set ESMF_DIR to the root of your ESMF installation.")
endif()

set(ESMF_MOD_DIR "${ESMF_DIR}/mod/modO/Linux.intel.64.intelmpi.default")
set(ESMF_LIB_DIR "${ESMF_DIR}/lib/libO/Linux.intel.64.intelmpi.default")
set(ESMF_INC_DIR "${ESMF_DIR}/include")

# Try to detect library name (adjust if your installation differs)
find_library(ESMF_LIBRARY
  NAMES esmf libesmf
  HINTS ${ESMF_LIB_DIR} ${ESMF_DIR}/lib
)

if (NOT ESMF_LIBRARY)
  message(FATAL_ERROR "Could not find ESMF library in ${ESMF_LIB_DIR}")
endif()

# -------------------------------------------------------------------------
# Build executable
add_executable(esmf_two_component_example
  ESMF_2component_fortran_example.f90
)

# Include directories
target_include_directories(esmf_two_component_example PRIVATE
  ${ESMF_INC_DIR}
  ${ESMF_MOD_DIR}
  ${MPI_Fortran_INCLUDE_PATH}
)

# Link libraries
target_link_libraries(esmf_two_component_example PRIVATE
  ${ESMF_LIBRARY}
  ${NETCDFF_LIBRARY} ${NETCDF_LIBRARY}
  MPI::MPI_Fortran
)

# Optionally show diagnostics
message(STATUS "NetCDF libraries: ${NETCDFF_LIBRARY} ${NETCDF_LIBRARY}")
message(STATUS "ESMF library: ${ESMF_LIBRARY}")
message(STATUS "ESMF include: ${ESMF_INC_DIR}")
message(STATUS "ESMF module:  ${ESMF_MOD_DIR}")
message(STATUS "MPI library:  ${MPI_Fortran_LIBRARIES}")

# -------------------------------------------------------------------------
# Build type and flags
set_target_properties(esmf_two_component_example PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules
)

# Example run instructions:
#   mkdir build && cd build
#   cmake -DESMF_DIR=/opt/esmf-8.6 .
#   make
#   mpirun -np 4 ./esmf_two_component_example

