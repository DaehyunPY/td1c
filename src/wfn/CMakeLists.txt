add_library(
    wfn STATIC
        clwfn.cpp
        wfn_orth.F90
        wfn_proj.F90
        wfn_projg.F90
        wfn_mask.F90
        wfn_ladapt.F90
        wfn_madapt.F90
        wfn_get_nradfc.F90
)
target_link_libraries(
    wfn PUBLIC
        OpenMP::OpenMP_Fortran
)