!######################################################################
subroutine hprod_kickk(knorm, wfn)

  use, intrinsic :: iso_c_binding
  use mod_const, only : iunit
  use mod_rad, only : xrad, nrad
  use mod_sph, only : nlat, cost
  use mod_ormas, only : nfcore, nfun

  implicit none
  real(c_double), intent(in) :: knorm
  complex(c_double_complex), intent(inout) :: wfn(1:*)

  integer(c_long) :: ifun
  integer(c_long) :: llr, ulr, irad, ilat
  complex(c_double_complex) :: faclv
  complex(c_double_complex), allocatable :: wfng(:,:,:)

  allocate(wfng(1:(nrad-1), 1:nlat, 1:nfun))
  call bas_sph2ang1_dyn(wfn, wfng)
  !$omp parallel default(shared) private(faclv, llr, ulr)
  !###########################
  call util_omp_disp(1, nrad-1, llr, ulr)
  do ilat = 1, nlat
     do irad = llr, ulr
        faclv = exp(iunit * knorm * xrad(irad) * cost(ilat))
        do ifun = nfcore + 1, nfun
           wfng(irad, ilat, ifun) = wfng(irad, ilat, ifun) * faclv
        end do
     end do
  end do
  !###########################
  !$omp end parallel
  call bas_ang2sph1_dyn(wfng, wfn)
  deallocate(wfng)

end subroutine hprod_kickk
!######################################################################
