!######################################################################
subroutine adi_v_cayley(icomp, dtime, lfield, wfn)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : ngrid, grid
  use mod_ormas, only : nfun

  implicit none
  integer(c_int), intent(in) :: icomp
  real(c_double), intent(in) :: dtime
  real(c_double), intent(in) :: lfield(1:3)
  complex(c_double_complex), intent(inout) :: wfn(1:ngrid, 1:nfun)

  integer(c_int) :: llg, ulg

  !$omp parallel default(shared) private(llg, ulg)
  !###########################
  call util_omp_disp(1, ngrid, llg, ulg)
  call adi_v_cayleyp(icomp, dtime, lfield, grid, wfn, llg, ulg)
  !###########################
  !$omp end parallel

end subroutine adi_v_cayley
!######################################################################
subroutine adi_v_cayleyp(icomp, dtime, lfield, grid, wfn, llg, ulg)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : ngrid
  use mod_ormas, only : nfcore, nfun
  use mod_const, only : one, two, runit, iunit

  implicit none
  integer(c_int), intent(in) :: icomp, llg, ulg
  real(c_double), intent(in) :: dtime
  real(c_double), intent(in) :: lfield(1:3)
  real(c_double), intent(in) :: grid(1:ngrid, 1:4)
  complex(c_double_complex), intent(inout) :: wfn(1:ngrid, 1:nfun)

  integer(c_int) :: ifun, igrid
  complex(c_double_complex) :: fac, ldip, vcay

  if (icomp == 1) then
     fac = iunit * dtime / two * lfield(3)
  else
     fac = runit * dtime / two * lfield(3)
  end if

  do ifun = nfcore + 1, nfun
     do igrid = llg, ulg
        ldip = fac * grid(igrid, 3)
        vcay = (one - ldip) / (one + ldip)
        wfn(igrid, ifun) = wfn(igrid, ifun) * vcay
     end do
  end do

end subroutine adi_v_cayleyp
!######################################################################
