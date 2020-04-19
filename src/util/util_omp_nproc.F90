!###########################################################
integer(c_int) function util_omp_nproc()

  use, intrinsic :: iso_c_binding

  implicit none
  ! local
  integer(c_int) :: iproc, nproc
  integer(c_int), external :: omp_get_num_threads
  integer(c_int), external :: omp_get_thread_num

  !$omp parallel default(shared) private(iproc)
  iproc = omp_get_thread_num()
  if (iproc == 0) nproc = omp_get_num_threads()
  !$omp end parallel

! DEBUG
! write(6, "('util_omp_nproc: nproc = ', i5)") nproc
! DEBUG

  util_omp_nproc = nproc
  return

end function util_omp_nproc
!###########################################################
