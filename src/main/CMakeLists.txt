add_executable(td1c td1c.cpp guess.cpp init.cpp tdse.cpp)
target_link_libraries(
    td1c PUBLIC
        # Targets in the project
        module io control mpi bas fedvr dvr sph tdcc ormas wfn
        hprod prop field gbas pes surff lapack util futil
        # External libraries
        rt pthread m dl  # System
        ifcore intlc  # Intel Compilers
        mkl_intel_lp64 mkl_core mkl_sequential
        gsl fgsl
        xcf03 xcf90 xc
        SHTOOLS-mp SHTOOLS
        ALGLIB
        OpenMP::OpenMP_C OpenMP::OpenMP_CXX OpenMP::OpenMP_Fortran
)