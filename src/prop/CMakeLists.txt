add_library(
    prop STATIC
        clh1prop.cpp
        clh2prop.cpp
        clcrk.cpp
        clirk.cpp
        clh1adi.cpp
        clh1exp.cpp
        clh1itr.cpp
        cletdrk.cpp
        cletdrb.cpp
        cletdrk2.cpp
        cletdrk4.cpp
        cllawrk.cpp
        cllawmul.cpp
        clcrnic.cpp
        cldpade.cpp
        clh1rat.cpp
        clexpab.cpp
        cllawab.cpp
        adi_gen_tadi.F90
        adi_t_explicit.F90
        adi_t_implicit.F90
        adi_laser_lgauge.F90
        adi_laser_vgauge.F90
        adi_v_cayley.F90
        h1exp_prop.F90
        h1exp_krylov.F90
        h1exp_exporb.F90
        h1itr_gen_cninv.F90
        h1itr_cninv.F90
        h1itr_prop.F90
        h1itr_prop2.F90
        dpade_gen.F90
        h1rat_gen.F90
        dpade_prod.F90
        h1rat_prod.F90
        h1rat_prodci.F90
        h1rat_prodci1.F90
        h1rat_prodci2.F90
        dpade_cnic.F90
)
target_link_libraries(
    prop PUBLIC
        hprod
        OpenMP::OpenMP_Fortran
)