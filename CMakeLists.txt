cmake_minimum_required(VERSION 3.12)
project(td1c LANGUAGES C CXX Fortran)

set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -stand=f18")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
list(APPEND CMAKE_INSTALL_RPATH "." "$ORIGIN" "$ORIGIN/lib" "$ORIGIN/../lib")
if(NOT "CMAKE_INSTALL_PREFIX" STREQUAL "")
    list(APPEND CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif()

# Ext packages
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

find_package(MKL REQUIRED)
include_directories("${MKL_ROOT}/include")
include_directories("${MKL_ROOT}/include/fftw")
link_directories("${MKL_ROOT}/lib/intel64")
list(APPEND CMAKE_INSTALL_RPATH "${MKL_ROOT}/lib/intel64")

find_package(GSL REQUIRED)
include_directories("${GSL_ROOT}/include")
link_directories("${GSL_ROOT}/lib")
list(APPEND CMAKE_INSTALL_RPATH "${GSL_ROOT}/lib")

find_package(FGSL REQUIRED)
include_directories("${FGSL_ROOT}/include/fgsl")
link_directories("${FGSL_ROOT}/lib")
list(APPEND CMAKE_INSTALL_RPATH "${FGSL_ROOT}/lib")

find_package(Libxc REQUIRED)
include_directories("${Libxc_ROOT}/include")
link_directories("${Libxc_ROOT}/lib")

find_package(SHTOOLS REQUIRED)
include_directories("${SHTOOLS_ROOT}/include")
link_directories("${SHTOOLS_ROOT}/lib")

find_package(ALGLIB REQUIRED)
include_directories("${ALGLIB_ROOT}/include")
link_directories("${ALGLIB_ROOT}/lib")

find_package(OpenMP REQUIRED)

list(APPEND TEMP_C "${CMAKE_C_IMPLICIT_LINK_DIRECTORIES}")
list(REMOVE_ITEM TEMP_C
    "${CMAKE_CXX_IMPLICIT_LINK_DIRECTORIES}"
    "${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES}")
list(APPEND TEMP_CXX "${CMAKE_CXX_IMPLICIT_LINK_DIRECTORIES}")
list(REMOVE_ITEM TEMP_CXX
    "${CMAKE_C_IMPLICIT_LINK_DIRECTORIES}"
    "${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES}")
list(APPEND TEMP_Fortran "${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES}")
list(REMOVE_ITEM TEMP_Fortran
    "${CMAKE_C_IMPLICIT_LINK_DIRECTORIES}"
    "${CMAKE_CXX_IMPLICIT_LINK_DIRECTORIES}")
list(APPEND CMAKE_INSTALL_RPATH
    "${TEMP_C}" "${TEMP_CXX}" "${TEMP_Fortran}"
    "${CMAKE_C_IMPLICIT_LINK_DIRECTORIES}"
    "${CMAKE_CXX_IMPLICIT_LINK_DIRECTORIES}"
    "${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES}")
list(REMOVE_DUPLICATES CMAKE_INSTALL_RPATH)

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