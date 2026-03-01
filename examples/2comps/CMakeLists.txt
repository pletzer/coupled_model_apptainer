# CMakeLists.txt — build the ESMF two-component coupling example
cmake_minimum_required(VERSION 3.18)
project(ESMF_TwoComponentExample LANGUAGES Fortran)

# --- User configuration section ------------------------------------------
# Set this to your ESMF installation path if not found automatically.
# For example:
#   cmake -DESMF_DIR=/path/to/esmf .
set(ESMF_LIBRARY "/software//lib/libO/Linux.intel.64.intelmpi.default/libesmf.a")
set(PIOC_LIBRARY "/software//lib/libO/Linux.intel.64.intelmpi.default/libpioc.a")

# NetCDF
set(NETCDF_LIBRARY "/software/lib/libnetcdf.a")
set(PNETCDF_LIBRARY "/software/lib/libpnetcdf.a")
set(NETCDFF_LIBRARY "/software/lib/libnetcdff.a")
set(HDF5_LIBRARY "/software/lib/libhdf5.a")
set(HDF5_HL_LIBRARY "/software/lib/libhdf5_hl.a") 


# MKL
#list(APPEND CMAKE_MODULE_PATH "$ENV{MKLROOT}/lib/cmake/mkl")
set(MKL_INTERFACE lp64)
find_package(MKL CONFIG REQUIRED PATHS $ENV{MKLROOT})

# -------------------------------------------------------------------------
# Find MPI (required by ESMF)
find_package(MPI REQUIRED Fortran)

set(ESMF_MOD_DIR "/software/mod/modO/Linux.intel.64.intelmpi.default")
set(ESMF_INC_DIR "/software/include")

# -------------------------------------------------------------------------
# Build esmfApp

# copy the namelist 
#configure_file("2comp_time_example.nml" "2comp_time_example.nml" COPYONLY) 
add_executable(esmfApp
  esmfApp.F90 esmf.F90 
  atm.F90 ocn.F90 conn.F90 vtk_mod.F90
)

# Include directories
target_include_directories(esmfApp PRIVATE
  ${ESMF_INC_DIR}
  ${ESMF_MOD_DIR}
  ${MPI_Fortran_INCLUDE_PATH}
)

# Link libraries
target_link_libraries(esmfApp PRIVATE
  ${ESMF_LIBRARY}
  ${PIOC_LIBRARY}
  ${PNETCDF_LIBRARY} ${NETCDFF_LIBRARY} ${NETCDF_LIBRARY}
  ${HDF5_HL_LIBRARY} ${HDF5_LIBRARY}
  MPI::MPI_Fortran
  MKL::MKL
  curl xml2 z stdc++
)
#==============================================================================

# Optionally show diagnostics
message(STATUS "NetCDF libraries: ${PNETCDF_LIBRARY} ${NETCDFF_LIBRARY} ${NETCDF_LIBRARY}")
message(STATUS "ESMF library: ${ESMF_LIBRARY}")
message(STATUS "PIOC library: ${PIOC_LIBRARY}")
message(STATUS "ESMF include: ${ESMF_INC_DIR}")
message(STATUS "ESMF module:  ${ESMF_MOD_DIR}")
message(STATUS "MPI library:  ${MPI_Fortran_LIBRARIES}")

