add_library(
    pes STATIC
        clpes.cpp
        pes_bind.F90
        pes_spec1_k.F90
        pes_spec1_kz.F90
)
target_link_libraries(
    pes PUBLIC
        OpenMP::OpenMP_C OpenMP::OpenMP_CXX OpenMP::OpenMP_Fortran
)