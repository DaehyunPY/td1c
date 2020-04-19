add_library(
    mpi STATIC
        clmpi.cpp
)
target_link_libraries(
    mpi PUBLIC
        OpenMP::OpenMP_C OpenMP::OpenMP_CXX
)