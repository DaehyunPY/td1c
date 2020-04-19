!######################################################################
subroutine hprod_poisson1(r2sph, v2sph)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad
  use mod_sph, only : lmax2

  implicit none
  complex(c_double_complex), intent(in) ::  r2sph(1:(nrad-1), 0:lmax2)
  complex(c_double_complex), intent(out) :: v2sph(1:(nrad-1), 0:lmax2)
  integer(c_long) :: lll, ull

  !$omp parallel default(shared) private(lll, ull)
  !###########################
  call util_omp_disp(0, lmax2, lll, ull)
  call hprod_poisson1p(r2sph, v2sph, lll, ull)
  !###########################
  !$omp end parallel

end subroutine hprod_poisson1
!######################################################################
subroutine hprod_poisson1p(rho2, v2sph, lll, ull)

  use, intrinsic :: iso_c_binding
  use mod_sph, only : lmax2
  use mod_rad, only : ndvr, nrad, xrad
  use mod_const, only : zero, one, czero, iunit
  use mod_bas, only : d2ll, bas_d2fac1, bas_d2fac2, bas_d2invr, bas_d2rpl0, bas_d2rpl1

  implicit none
  integer(c_long), intent(in) :: lll, ull
  complex(c_double_complex), intent(in) ::     rho2(1:(nrad-1), 0:lmax2)
  complex(c_double_complex), intent(inout) :: v2sph(1:(nrad-1), 0:lmax2)

  real(c_double) :: tmpr, tmpi
  real(c_double), allocatable :: rrho2(:,:)
  integer(c_long) :: dim, l, irad, ld, info

  ld = ndvr + 1
  dim = nrad - 1
  allocate(rrho2(1:dim, 1:2))
  do l = 0, lmax2
     bas_d2fac1(l) = one / xrad(dim + 1) ** (2 * l + 1)
  end do

  do l = lll, ull
     ! r.h.s and BC
     tmpr = zero
     tmpi = zero
     do irad = 1, dim
        rrho2(irad, 1) = dble( rho2(irad, l))
        rrho2(irad, 2) = aimag(rho2(irad, l))
        tmpr = tmpr + rrho2(irad, 1) * bas_d2rpl0(irad, l)
        tmpi = tmpi + rrho2(irad, 2) * bas_d2rpl0(irad, l)
        rrho2(irad, 1) = rrho2(irad, 1) * bas_d2invr(irad, l)
        rrho2(irad, 2) = rrho2(irad, 2) * bas_d2invr(irad, l)
     end do
     tmpr = tmpr * bas_d2fac1(l)
     tmpi = tmpi * bas_d2fac1(l)
     ! Poisson's equation with vanishing BCs
     call DPBTRS('L', dim, ndvr, 2, d2ll(1,1,l), ld, rrho2, dim, info)
     if (info /= 0) then
        write(6, "('hprod_poisson1p: dpbtrs bad info. ', 2i5)") l, info
        stop
     end if
     ! solution with the correct BC
     do irad = 1, dim
        v2sph(irad, l) = (rrho2(irad, 1) + tmpr * bas_d2rpl1(irad, l)) * bas_d2fac2(irad, l) &
                       + (rrho2(irad, 2) + tmpi * bas_d2rpl1(irad, l)) * bas_d2fac2(irad, l) * iunit
     end do
  end do

  deallocate(rrho2)

end subroutine hprod_poisson1p
!######################################################################
