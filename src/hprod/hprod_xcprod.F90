!######################################################################
subroutine hprod_xcprod(dofc, wfn, hwfn)

  ! dynamical-core-eXchange

  use, intrinsic :: iso_c_binding
  use mod_sph, only : nlat
  use mod_ormas, only : nfcore, nocc, nfun
  use mod_rad, only : nrad, nradfc, xrad, wrad
  use mod_const, only : zero, two, ctwo, pi

  implicit none
  logical(c_bool), intent(in) :: dofc
  complex(c_double_complex), intent(inout) :: wfn(1:(nrad-1), 1:nlat, 1:nfun)
  complex(c_double_complex), intent(inout) :: hwfn(1:(nrad-1), 1:nlat, 1:nfun)

  integer(c_int) :: lll, ull, iproc
  integer(c_int) :: ifun, ilat, irad
  real(c_double), allocatable :: rho1(:,:)
  real(c_double), allocatable :: v2xc(:,:)
  integer(c_int), external :: util_omp_iproc
  
  allocate(rho1(1:(nrad-1), 1:nlat))
  allocate(v2xc(1:(nrad-1), 1:nlat))
  rho1(1:(nrad-1), 1:nlat) = zero
  v2xc(1:(nrad-1), 1:nlat) = zero

  !$omp parallel default(shared) private(iproc, lll, ull)
  !###########################
  iproc = util_omp_iproc()
  call util_omp_disp(1, nlat, lll, ull)
  do ifun = 1, nocc
     do ilat = lll, ull
        do irad = 1, nrad - 1
           rho1(irad, ilat) = rho1(irad, ilat) + conjg(wfn(irad, ilat, ifun)) &
                                                     * wfn(irad, ilat, ifun)
        end do
     end do
  end do
  do ilat = lll, ull
     do irad = 1, nrad - 1
!       rho1(irad, ilat) =  ctwo * rho1(irad, ilat) / (xrad(irad) * xrad(irad) * wrad(irad))
!       rho1(irad, ilat) =  ctwo * rho1(irad, ilat) / (xrad(irad) * xrad(irad) * wrad(irad)) / (two * pi)
        rho1(irad, ilat) =  ctwo * rho1(irad, ilat) / (xrad(irad) * xrad(irad)) / (two * pi)
     end do
  end do

  ! acting on FC
  if (dofc) then
     do ifun = 1, nfcore
        do ilat = lll, ull
           do irad = 1, nradfc
              hwfn(irad, ilat, ifun) = hwfn(irad, ilat, ifun) + v2xc(irad, ilat) * wfn(irad, ilat, ifun)
           end do
        end do
     end do
  end if
  ! acting on DC and active
  do ifun = nfcore + 1, nfun
     do ilat = lll, ull
        do irad = 1, nrad - 1
           hwfn(irad, ilat, ifun) = hwfn(irad, ilat, ifun) + v2xc(irad, ilat) * wfn(irad, ilat, ifun)
        end do
     end do
  end do
  !###########################
  !$omp end parallel

  deallocate(v2xc)
  deallocate(rho1)
  
end subroutine hprod_xcprod
!######################################################################
