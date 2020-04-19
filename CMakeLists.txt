cmake_minimum_required(VERSION 3.12)
project(td1c LANGUAGES C CXX Fortran)

set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -stand=f18")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
set(CMAKE_INSTALL_RPATH ".:$ORIGIN:$ORIGIN/lib:$ORIGIN/../lib")
if(NOT "CMAKE_INSTALL_PREFIX" STREQUAL "")
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_RPATH}:${CMAKE_INSTALL_PREFIX}/lib")
endif()
if(NOT "$ENV{LD_LIBRARY_PATH}" STREQUAL "")
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_RPATH}:$ENV{LD_LIBRARY_PATH}")
endif()

# Ext packages
find_package(OpenMP REQUIRED)  # https://www.openmp.org/

# Ext package - MKL https://software.intel.com/mkl/
if(NOT "$ENV{MKLROOT}" STREQUAL "")
    message("-- Env MKLROOT: $ENV{MKLROOT}")
    include_directories("$ENV{MKLROOT}/include")
    include_directories("$ENV{MKLROOT}/include/fftw")
    link_directories("$ENV{MKLROOT}/lib/intel64")
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_RPATH}:$ENV{MKLROOT}/lib/intel64")
else()
    message("-- Env MKLROOT: unset")
endif()

# Ext package - ALGLIB https://www.alglib.net/
include_directories(/opt/lab/alglib/inc)
link_directories(/opt/lab/alglib/obj)

# Ext package - GSL https://www.gnu.org/software/gsl/
include_directories(/opt/lab/gsl/include)
include_directories(/opt/lab/gsl/include/fgsl)
link_directories(/opt/lab/gsl/lib)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_RPATH}:/opt/lab/gsl/lib")

# Ext package - Libxc https://www.tddft.org/programs/libxc/
include_directories(/opt/lab/libxc/include)
link_directories(/opt/lab/libxc/lib)

# Ext package - SHTOOLS https://shtools.oca.eu/shtools/public/
include_directories(/opt/lab/SHTOOLS/modules)
link_directories(/opt/lab/SHTOOLS/lib)

# Libraries
set(CMAKE_Fortran_MODULE_DIRECTORY "${PROJECT_SOURCE_DIR}/include")
include_directories("${PROJECT_SOURCE_DIR}/include")
add_subdirectory(src/module)
add_subdirectory(src/io)
add_subdirectory(src/control)
add_subdirectory(src/mpi)
add_subdirectory(src/bas)
add_subdirectory(src/fedvr)
add_subdirectory(src/dvr)
add_subdirectory(src/sph)
add_subdirectory(src/tdcc)
add_subdirectory(src/ormas)
add_subdirectory(src/wfn)
add_subdirectory(src/hprod)
add_subdirectory(src/prop)
add_subdirectory(src/field)
add_subdirectory(src/gbas)
add_subdirectory(src/pes)
add_subdirectory(src/surff)
add_subdirectory(src/lapack)
add_subdirectory(src/util)
add_subdirectory(src/futil)

# Executables
add_subdirectory(src/main)

install(
    TARGETS module bas tdcc ormas td1c
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
)