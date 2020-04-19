add_library(
    lapack STATIC
        lapack_sv.F90
        lapack_ev.F90
        zclear_omp.F90
        zcopy_omp.F90
        zscal_omp.F90
        zxpy_omp.F90
        zxpyz_omp.F90
        zxmyz_omp.F90
        zaxpy_omp.F90
        zaxpyz_omp.F90
        zaxpbyz_omp.F90
        zdotc_omp.F90
        zadd_omp.F90
)
target_link_libraries(
    lapack PUBLIC
        OpenMP::OpenMP_Fortran
)