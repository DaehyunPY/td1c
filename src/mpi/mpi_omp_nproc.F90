!###########################################################
integer(c_int) function mpi_omp_nproc()

  use, intrinsic :: iso_c_binding

  implicit none
  integer(c_int) :: iproc, nproc
  integer(c_int), external :: omp_get_num_threads
  integer(c_int), external :: util_omp_iproc

  !$omp parallel default(shared) private(iproc)
  iproc = util_omp_iproc()
  if (iproc == 0) nproc = omp_get_num_threads()
  !$omp end parallel

! DEBUG
! write(6, "('mpi_omp_nproc: nproc = ', i5)") nproc
! DEBUG

  mpi_omp_nproc = nproc
  return

end function mpi_omp_nproc
!###########################################################
