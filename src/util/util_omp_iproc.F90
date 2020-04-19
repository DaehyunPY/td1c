!###########################################################
integer(c_long) function util_omp_iproc()

  use, intrinsic :: iso_c_binding

  implicit none
  ! local
  integer, external :: omp_get_thread_num
  util_omp_iproc = omp_get_thread_num()

end function util_omp_iproc
!###########################################################
