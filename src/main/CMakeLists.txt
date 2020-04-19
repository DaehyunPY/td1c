add_executable(td1c td1c.cpp guess.cpp init.cpp tdse.cpp)
target_link_libraries(
    td1c PUBLIC
        # Targets in the project
        module io control mpi bas fedvr dvr sph tdcc ormas wfn
        hprod prop field gbas pes surff lapack util futil
        # External libraries
        SHTOOLS
        fgsl gsl
        xcf03 xcf90 xc
        mkl_intel_lp64 mkl_core mkl_sequential ifcore
        rt pthread m dl
        /opt/lab/alglib/obj/alglibinternal.o
        /opt/lab/alglib/obj/alglibmisc.o
        /opt/lab/alglib/obj/ap.o
        /opt/lab/alglib/obj/dataanalysis.o
        /opt/lab/alglib/obj/diffequations.o
        /opt/lab/alglib/obj/fasttransforms.o
        /opt/lab/alglib/obj/integration.o
        /opt/lab/alglib/obj/interpolation.o
        /opt/lab/alglib/obj/linalg.o
        /opt/lab/alglib/obj/optimization.o
        /opt/lab/alglib/obj/solvers.o
        /opt/lab/alglib/obj/specialfunctions.o
        /opt/lab/alglib/obj/statistics.o
)