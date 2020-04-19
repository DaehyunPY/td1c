!######################################################################
! test FC use nfcore_tdcis instead of nfcore
! test multipole expansion
!######################################################################
subroutine hprod_mkv2mf_tdcis_init()

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nfun, nfcore1, nfcore2, nfcore
  use mod_rad, only : nrad, nradfc, nradgs
  use mod_sph, only : lmax2
  use mod_hprod, only : rho2, v2sph, v2ang

  implicit none
  integer(c_int) :: lll, ull

  !$omp parallel default(shared) private(lll, ull)
  !###########################
  call util_omp_disp(0, lmax2, lll, ull)
  ! original full Poisson [0, nrad]
  ! call hprod_mkv2mf_poissonp_tdcis(rho2, v2sph, lll, ull)
  ! Poisson [0, nradgs] and multipole expansion (nradgs, nrad]
  call hprod_mkv2mf_poissonp_tdcis_init(rho2, v2sph, lll, ull)
  !###########################
  !$omp end parallel

  call bas_sph2ang2(v2sph, v2ang)

end subroutine hprod_mkv2mf_tdcis_init
!######################################################################
subroutine hprod_mkv2mf_tdcis()

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nfun, nfcore1, nfcore2, nfcore
  use mod_rad, only : nrad, nradfc, nradgs
  use mod_sph, only : lmax2
  use mod_hprod, only : rho2, v2sph, v2ang

  implicit none
  integer(c_int) :: lll, ull

  !$omp parallel default(shared) private(lll, ull)
  !###########################
  call util_omp_disp(0, lmax2, lll, ull)
  call hprod_mkv2mf_poissonp_tdcis(rho2, v2sph, lll, ull)
  !###########################
  !$omp end parallel

  call hprod_mkv2mf_herm_tdcis(v2sph)
  call bas_sph2ang2(v2sph, v2ang)

end subroutine hprod_mkv2mf_tdcis
!######################################################################
subroutine hprod_mkv2mf_poissonp_tdcis(rho2, v2sph, lll, ull)

  use, intrinsic :: iso_c_binding
  use mod_sph, only : lmax2
  use mod_ormas, only : nfun, nfcore
  use mod_rad, only : ndvr, nrad, xrad
  use mod_const, only : zero, one, czero, iunit
  use mod_bas, only : mval, d2ll, bas_d2fac1, bas_d2fac2, bas_d2invr, bas_d2rpl0, bas_d2rpl1

  use mod_ormas, only : nfcore_tdcis
  use mod_rad, only : nradgs, ecs_flag, irad_ecs
  use mod_const, only : one, half, four, PI

  implicit none
  integer(c_int), intent(in) :: lll, ull
  complex(c_double_complex), intent(in) ::     rho2(1:(nrad-1), 0:lmax2, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(inout) :: v2sph(1:(nrad-1), 0:lmax2, 1:nfun, 1:nfun)

  real(c_double) :: tmpr, tmpi
  real(c_double), allocatable :: rrho2(:,:)
  integer(c_int) :: ifun, jfun, l, m, irad, ld, info
  integer(c_int) :: dim, dimbc

  ! multipole expansion to be removed
  ! complex(c_double_complex) :: charge
  ! complex(c_double_complex) :: mpfac
  ! multipole expansion to be removed

  ld = ndvr + 1

  dim = nrad - 1
  dimbc = min(nradgs, nrad - 1)
! tdcis_sato
  if (ecs_flag == 1) dimbc = min(dimbc, irad_ecs-1)
! tdcis_sato

  allocate(rrho2(1:dim, 1:2))

  do l = lll, ull
     bas_d2fac1(l) = one / xrad(dimbc + 1) ** (2 * l + 1)
  end do

  do ifun = nfcore_tdcis + 1, nfun
     do jfun = ifun, nfun
        !  do ifun = 1, nfun
        !     do jfun = 1, nfun
        do l = lll, ull
           ! multipole expansion
           ! mpfac = four*PI/(2 * l + 1) 
           ! multipole expansion
           m = -mval(ifun) + mval(jfun)
           if (l < abs(m)) then
              v2sph(1:dim, l, jfun, ifun) = czero 
              cycle
           end if
           ! r.h.s and BC
           tmpr = zero
           tmpi = zero
           do irad = 1, dimbc
              rrho2(irad, 1) = dble( rho2(irad, l, jfun, ifun))
              rrho2(irad, 2) = aimag(rho2(irad, l, jfun, ifun))
              tmpr = tmpr + rrho2(irad, 1) * bas_d2rpl0(irad, l)
              tmpi = tmpi + rrho2(irad, 2) * bas_d2rpl0(irad, l)
              rrho2(irad, 1) = rrho2(irad, 1) * bas_d2invr(irad, l)
              rrho2(irad, 2) = rrho2(irad, 2) * bas_d2invr(irad, l)
           end do
           tmpr = tmpr * bas_d2fac1(l)
           tmpi = tmpi * bas_d2fac1(l)
           ! Poisson's equation with vanishing BCs
           call DPBTRS('L', dimbc, ndvr, 2, d2ll(1,1,l), ld, rrho2, dim, info)
           if (info /= 0) then
              write(6, "('hprod_poisson: dpbtrs bad info. ', 5i5)") ifun, jfun, l, m, info
              stop
           end if
           ! solution with the correct BC
           do irad = 1, dimbc
              v2sph(irad, l, jfun, ifun) = (rrho2(irad, 1) + tmpr * bas_d2rpl1(irad, l)) * bas_d2fac2(irad, l) &
                                         + (rrho2(irad, 2) + tmpi * bas_d2rpl1(irad, l)) * bas_d2fac2(irad, l) * iunit
           end do
           ! multipole expansion -> to be removed
           !    charge = (tmpr + tmpi * iunit)/bas_d2fac1(l)
           !    do irad = dimbc + 1, dim
           !       v2sph(irad, l, jfun, ifun) = mpfac / xrad(irad)**(l + 1) * charge
           !    end do
        end do
     end do
  end do

  deallocate(rrho2)

end subroutine hprod_mkv2mf_poissonp_tdcis
!######################################################################
subroutine hprod_mkv2mf_poissonp_tdcis_init(rho2, v2sph, lll, ull)

  use, intrinsic :: iso_c_binding
  use mod_sph, only : lmax2
  use mod_ormas, only : nfun, nfcore
  use mod_rad, only : ndvr, nrad, xrad
  use mod_const, only : zero, one, czero, iunit
  use mod_bas, only : mval, d2ll, bas_d2fac1, bas_d2fac2, bas_d2invr, bas_d2rpl0, bas_d2rpl1

  use mod_ormas, only : nfcore_tdcis
  use mod_rad, only : nradgs, ecs_flag, irad_ecs
  use mod_const, only : one, half, four, PI

  implicit none
  integer(c_int), intent(in) :: lll, ull
  complex(c_double_complex), intent(in) ::     rho2(1:(nrad-1), 0:lmax2, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(inout) :: v2sph(1:(nrad-1), 0:lmax2, 1:nfun, 1:nfun)

  real(c_double) :: tmpr, tmpi
  real(c_double), allocatable :: rrho2(:,:)
  integer(c_int) :: ifun, jfun, l, m, irad, ld, info
  integer(c_int) :: dim, dimbc

  ! multipole expansion
  complex(c_double_complex) :: charge
  complex(c_double_complex) :: mpfac
  ! multipole expansion

  ld = ndvr + 1

  dim = nrad - 1
  dimbc = min(nradgs, nrad - 1)
! tdcis_sato
  if (ecs_flag == 1) dimbc = min(dimbc, irad_ecs-1)
! tdcis_sato

  allocate(rrho2(1:dim, 1:2))

  do l = lll, ull
     bas_d2fac1(l) = one / xrad(dimbc + 1) ** (2 * l + 1)
  end do

  do ifun = 1, nfun
     do jfun = 1, nfun
        do l = lll, ull
           ! multipole expansion
           mpfac = four*PI/(2 * l + 1) 
           ! multipole expansion
           m = -mval(ifun) + mval(jfun)
           if (l < abs(m)) then
              v2sph(1:dim, l, jfun, ifun) = czero 
              cycle
           end if
           ! r.h.s and BC
           tmpr = zero
           tmpi = zero
           do irad = 1, dimbc
              rrho2(irad, 1) = dble( rho2(irad, l, jfun, ifun))
              rrho2(irad, 2) = aimag(rho2(irad, l, jfun, ifun))
              tmpr = tmpr + rrho2(irad, 1) * bas_d2rpl0(irad, l)
              tmpi = tmpi + rrho2(irad, 2) * bas_d2rpl0(irad, l)
              rrho2(irad, 1) = rrho2(irad, 1) * bas_d2invr(irad, l)
              rrho2(irad, 2) = rrho2(irad, 2) * bas_d2invr(irad, l)
           end do
           tmpr = tmpr * bas_d2fac1(l)
           tmpi = tmpi * bas_d2fac1(l)
           ! Poisson's equation with vanishing BCs
           call DPBTRS('L', dimbc, ndvr, 2, d2ll(1,1,l), ld, rrho2, dim, info)
           if (info /= 0) then
              write(6, "('hprod_poisson: dpbtrs bad info. ', 5i5)") ifun, jfun, l, m, info
              stop
           end if
           ! solution with the correct BC
           do irad = 1, dimbc
              v2sph(irad, l, jfun, ifun) = (rrho2(irad, 1) + tmpr * bas_d2rpl1(irad, l)) * bas_d2fac2(irad, l) &
                                         + (rrho2(irad, 2) + tmpi * bas_d2rpl1(irad, l)) * bas_d2fac2(irad, l) * iunit
           end do
           ! multipole expansion
           charge = (tmpr + tmpi * iunit)/bas_d2fac1(l)
           do irad = dimbc + 1, dim
              v2sph(irad, l, jfun, ifun) = mpfac / xrad(irad)**(l + 1) * charge
           end do
        end do
     end do
  end do

  deallocate(rrho2)

end subroutine hprod_mkv2mf_poissonp_tdcis_init
!######################################################################
subroutine hprod_mkv2mf_herm_tdcis(v2sph)
  use, intrinsic :: iso_c_binding
  use mod_sph, only : lmax2, nlat
  use mod_rad, only : nrad, nradfc
  use mod_ormas, only : nfcore, nfcore1, nfcore2, nfun
  use mod_ormas, only : nfcore_tdcis
  use mod_rad, only : ecs_flag, irad_ecs

  implicit none 
  complex(c_double_complex), intent(inout) :: v2sph(1:(nrad-1), 0:lmax2, 1:nfun, 1:nfun)
  integer(c_int) :: l, ilat, ifun, jfun, irad
  integer(c_int) :: llr, ulr, nradmax

  if (ecs_flag == 0) then
     nradmax = nrad - 1
  else 
     nradmax = irad_ecs - 1
  end if

  !$omp parallel default(shared) private(llr, ulr)
  call util_omp_disp(1, nradmax, llr, ulr)
  do ifun = nfcore_tdcis + 1, nfun
     do jfun = ifun + 1, nfun
        do l = 0, lmax2
           do irad = llr, ulr
              v2sph(irad, l, ifun, jfun) = conjg(v2sph(irad, l, jfun, ifun))
           end do
        end do
     end do
  end do
  !$omp end parallel
end subroutine hprod_mkv2mf_herm_tdcis
!######################################################################
! subroutine hprod_mkv2mf_poissonp_tdcis(rho2, v2sph, lll, ull)

!   use, intrinsic :: iso_c_binding
!   use mod_sph, only : lmax2
!   use mod_ormas, only : nfun, nfcore
!   use mod_rad, only : ndvr, nrad, xrad
!   use mod_const, only : zero, one, czero, iunit
!   use mod_bas, only : mval, d2ll, bas_d2fac1, bas_d2fac2, bas_d2invr, bas_d2rpl0, bas_d2rpl1

!   use mod_ormas, only : nfcore_tdcis

!   implicit none
!   integer(c_int), intent(in) :: lll, ull
!   complex(c_double_complex), intent(in) ::     rho2(1:(nrad-1), 0:lmax2, 1:nfun, 1:nfun)
!   complex(c_double_complex), intent(inout) :: v2sph(1:(nrad-1), 0:lmax2, 1:nfun, 1:nfun)

!   real(c_double) :: tmpr, tmpi
!   real(c_double), allocatable :: rrho2(:,:)
!   integer(c_int) :: ifun, jfun, l, m, irad, ld, info
!   integer(c_int) :: dim
!   ld = ndvr + 1
!   dim = nrad - 1
!   allocate(rrho2(1:dim, 1:2))

!   do l = lll, ull
!      bas_d2fac1(l) = one / xrad(dim + 1) ** (2 * l + 1)
!   end do

!   ! do ifun = nfcore_tdcis + 1, nfun
!   !    do jfun = nfcore_tdcis + 1, nfun
!   do ifun = 1, nfun
!      do jfun = 1, nfun
!         !     do jfun = jfun_ll, jfun_ul
!         do l = lll, ull
!            m = -mval(ifun) + mval(jfun)
!            if (l < abs(m)) then
!               v2sph(1:dim, l, jfun, ifun) = czero 
!               cycle
!            end if
!            ! r.h.s and BC
!            tmpr = zero
!            tmpi = zero
!            do irad = 1, dim
!               rrho2(irad, 1) = dble( rho2(irad, l, jfun, ifun))
!               rrho2(irad, 2) = aimag(rho2(irad, l, jfun, ifun))
!               tmpr = tmpr + rrho2(irad, 1) * bas_d2rpl0(irad, l)
!               tmpi = tmpi + rrho2(irad, 2) * bas_d2rpl0(irad, l)
!               rrho2(irad, 1) = rrho2(irad, 1) * bas_d2invr(irad, l)
!               rrho2(irad, 2) = rrho2(irad, 2) * bas_d2invr(irad, l)
!            end do
!            tmpr = tmpr * bas_d2fac1(l)
!            tmpi = tmpi * bas_d2fac1(l)
!            ! Poisson's equation with vanishing BCs
!            call DPBTRS('L', dim, ndvr, 2, d2ll(1,1,l), ld, rrho2, dim, info)
!            if (info /= 0) then
!               write(6, "('hprod_poisson: dpbtrs bad info. ', 5i5)") ifun, jfun, l, m, info
!               stop
!            end if
!            ! solution with the correct BC
!            do irad = 1, dim
!               v2sph(irad, l, jfun, ifun) = (rrho2(irad, 1) + tmpr * bas_d2rpl1(irad, l)) * bas_d2fac2(irad, l) &
!                                          + (rrho2(irad, 2) + tmpi * bas_d2rpl1(irad, l)) * bas_d2fac2(irad, l) * iunit
!            end do
!         end do
!      end do
!   end do

!   deallocate(rrho2)

! end subroutine hprod_mkv2mf_poissonp_tdcis
!######################################################################
