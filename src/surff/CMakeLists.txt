add_library(
    surff STATIC
        integ_trap.cpp
        integ_simpson.cpp
        integ_CN.cpp
        integ_EI.cpp
        surff_init.cpp
        surff_prop.cpp
        surff_phys.cpp
        calc_orb.F90
        surff_rec_v2xmat.F90
)
target_link_libraries(
    surff PUBLIC
        OpenMP::OpenMP_C OpenMP::OpenMP_CXX OpenMP::OpenMP_Fortran
)