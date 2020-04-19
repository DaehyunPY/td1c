!###########################################################
subroutine mpi_omp_divide(ll, ul, llp, ulp)

  use, intrinsic :: iso_c_binding

  implicit none
  integer(c_int), intent(in) :: ll, ul
  integer(c_int), intent(out) :: llp, ulp

  integer(c_int) :: ndata, ndpp, iproc, nproc
  integer(c_int), external :: util_omp_nproc
  integer(c_int), external :: util_omp_iproc

  iproc = util_omp_iproc()
  nproc = util_omp_iproc()
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

end subroutine mpi_omp_divide
!###########################################################
