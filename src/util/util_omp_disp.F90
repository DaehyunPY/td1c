!###########################################################
subroutine util_omp_disp(ll, ul, llp, ulp)

  use, intrinsic :: iso_c_binding

  implicit none
  integer(c_long), intent(in) :: ll, ul
  integer(c_long), intent(out) :: llp, ulp
  ! local
  integer(c_long) :: ndata, ndpp, iproc, nproc
  integer, external :: omp_get_num_threads
  integer, external :: omp_get_thread_num

  iproc = omp_get_thread_num()
  nproc = omp_get_num_threads()
! DEBUG
!  write(6, "('iproc = ', i5)") iproc
!  write(6, "('nproc = ', i5)") nproc
! DEBUG

  ndata = ul - ll + 1
  ndpp = ndata / nproc
! DEBUG
!  write(6, "('ndata = ', i5)") ndata
!  write(6, "('ndpp  = ', i5)") ndpp
! DEBUG
  if (iproc == nproc - 1) then
     llp = ndpp * iproc + ll
     ulp = ul
  else
     llp = ndpp * iproc + ll
     ulp = llp + ndpp - 1
  end if

! DEBUG
!  write(6, "('iproc = ', i5, 'llp   = ', i5)") iproc, llp
!  write(6, "('iproc = ', i5, 'ulp   = ', i5)") iproc, ulp
! DEBUG

end subroutine util_omp_disp
!###########################################################
subroutine util_omp_disp_top(iproc, nproc, ll, ul, llp, ulp)

  use, intrinsic :: iso_c_binding

  implicit none
  integer(c_long), intent(in) :: iproc, nproc, ll, ul
  integer(c_long), intent(out) :: llp, ulp
  ! local
  integer(c_long) :: ndata, ndpp

  ndata = ul - ll + 1
  ndpp = ndata / nproc
! DEBUG
!  write(6, "('ndata = ', i5)") ndata
!  write(6, "('ndpp  = ', i5)") ndpp
! DEBUG
  if (iproc == nproc - 1) then
     llp = ndpp * iproc + ll
     ulp = ul
  else
     llp = ndpp * iproc + ll
     ulp = llp + ndpp - 1
  end if

! DEBUG
!  write(6, "('iproc = ', i5, 'llp   = ', i5)") iproc, llp
!  write(6, "('iproc = ', i5, 'ulp   = ', i5)") iproc, ulp
! DEBUG

end subroutine util_omp_disp_top
!###########################################################
