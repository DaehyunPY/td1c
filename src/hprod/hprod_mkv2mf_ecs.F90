!######################################################################
subroutine hprod_mkv2mf_ecs()

  use, intrinsic :: iso_c_binding
  use mod_hprod, only : rho2, v2sph, v2ang

  implicit none
  call hprod_mkv2mf_poisson_ecs('tot', rho2, v2sph)
  call bas_sph2ang2_ecs(v2sph, v2ang);
  !!call hprod_mkv2mf_herm_ecs('tot', v2sph, v2ang)

end subroutine hprod_mkv2mf_ecs
!######################################################################
subroutine hprod_mkv2mf_dyn_ecs()

  use, intrinsic :: iso_c_binding
  use mod_hprod, only : rho2, v2sph, v2ang
  use mod_hprod, only : rho23j, v2sph3j, v2ange, v2ango

  implicit none
  call hprod_mkv2mf_poisson_ecs('dyn', rho2, v2sph)
! Sato_ECS
  call bas_sph2ang2_dyn_ecs(v2sph, v2ang);
  !!call hprod_mkv2mf_herm_ecs('dyn', v2sph, v2ang)
!  call bas_sph2ang2_dyn(v2sph, v2ang);
!  call hprod_mkv2mf_herm('dyn', v2sph, v2ang)
! Sato_ECS

  call hprod_mkv2mf_poisson_ecs('fc1dyn', rho2, v2sph)
! Sato_ECS
  call bas_sph2ang2_fc1dyn_ecs(v2sph, v2ang);
  !!call hprod_mkv2mf_herm_ecs('fc1dyn', v2sph, v2ang)
!  call bas_sph2ang2_fc1dyn(v2sph, v2ang);
!  call hprod_mkv2mf_herm('fc1dyn', v2sph, v2ang)
! Sato_ECS
  
!def  call hprod_mkv2mf_poisson('dyn', rho2, v2sph)
!def  call bas_sph2ang2_dyn(v2sph, v2ang);
!def  call hprod_mkv2mf_herm('dyn', v2sph, v2ang)
!def  call hprod_mkv2mf_poisson('fc1dyn', rho2, v2sph)
!def  call bas_sph2ang2_fc1dyn(v2sph, v2ang);
!def  call hprod_mkv2mf_herm('fc1dyn', v2sph, v2ang)

end subroutine hprod_mkv2mf_dyn_ecs
!######################################################################
subroutine hprod_mkv2mf_poisson_ecs(v2_type, rho2, v2sph)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nfun, nfcore1, nfcore2, nfcore
  use mod_rad, only : nrad, nradfc
  use mod_sph, only : lmax2

  implicit none
  character(len=*), intent(in) :: v2_type
  complex(c_double_complex), intent(in) ::     rho2(1:(nrad-1), 0:lmax2, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(inout) :: v2sph(1:(nrad-1), 0:lmax2, 1:nfun, 1:nfun)
  integer(c_long) :: dim, ifun_ll, ifun_ul, jfun_ll, jfun_ul, lll, ull

  if (trim(v2_type) == 'tot') then
     dim = nrad - 1
     ifun_ll = 1
     ifun_ul = nfun
     jfun_ll = 1
     jfun_ul = nfun
  else if (trim(v2_type) == 'fcfc') then
!BUG dim = nradfc
     dim = nrad - 1
     ifun_ll = 1
     ifun_ul = nfcore
     jfun_ll = 1
     jfun_ul = nfcore
  else if (trim(v2_type) == 'fc1fc1') then
!BUG dim = nradfc
     dim = nrad - 1
     ifun_ll = nfcore2 + 1
     ifun_ul = nfcore
     jfun_ll = nfcore2 + 1
     jfun_ul = nfcore
  else if (trim(v2_type) == 'fc2fc2') then
!BUG dim = nradfc
     dim = nrad - 1
     ifun_ll = 1
     ifun_ul = nfcore2
     jfun_ll = 1
     jfun_ul = nfcore2
  else if (trim(v2_type) == 'fc1dyn') then
     dim = nradfc
     ifun_ll = nfcore2 + 1
     ifun_ul = nfcore
     jfun_ll = nfcore + 1
     jfun_ul = nfun
  else if (trim(v2_type) == 'fc2dyn') then
     dim = nradfc
     ifun_ll = 1
     ifun_ul = nfcore2
     jfun_ll = nfcore + 1
     jfun_ul = nfun
  else if (trim(v2_type) == 'dyn') then
     dim = nrad - 1
     ifun_ll = nfcore + 1
     ifun_ul = nfun
     jfun_ll = nfcore + 1
     jfun_ul = nfun
  else
     stop "bad v2_type in hprod_mkv2mf_poisson"
  end if

  !$omp parallel default(shared) private(lll, ull)
  !###########################
  call util_omp_disp(0, lmax2, lll, ull)
  call hprod_mkv2mf_poissonp_ecs(rho2, v2sph, dim, lll, ull, ifun_ll, ifun_ul, jfun_ll, jfun_ul)
  !------------------------------------------------------------
  !NEW orimo
  if (trim(v2_type) == 'fc1dyn' .or. trim(v2_type) == 'fc2dyn') then
     call hprod_mkv2mf_poissonp_ecs(rho2, v2sph, dim, lll, ull, jfun_ll, jfun_ul, ifun_ll, ifun_ul)
  end if
  !------------------------------------------------------------
  !###########################
  !$omp end parallel


end subroutine hprod_mkv2mf_poisson_ecs
!######################################################################
subroutine hprod_mkv2mf_poissonp_ecs(rho2, v2sph, dim, lll, ull, ifun_ll, ifun_ul, jfun_ll, jfun_ul)

  use, intrinsic :: iso_c_binding
  use mod_sph, only : lmax2
  use mod_ormas, only : nfun
  use mod_rad, only : ndvr, nrad, xrad, irad_ecs
  use mod_const, only : zero, one, czero, iunit
  use mod_bas, only : mval, d2ll, bas_d2fac1, bas_d2fac2, bas_d2invr, bas_d2rpl0, bas_d2rpl1, bas_d2crpl1

  implicit none
  integer(c_long), intent(in) :: dim, lll, ull, ifun_ll, ifun_ul, jfun_ll, jfun_ul
  complex(c_double_complex), intent(in) ::     rho2(1:(nrad-1), 0:lmax2, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(inout) :: v2sph(1:(nrad-1), 0:lmax2, 1:nfun, 1:nfun)

  real(c_double) :: tmpr, tmpi
  real(c_double), allocatable :: rrho2(:,:)
  integer(c_long) :: ifun, jfun, l, m, irad, ld, info

  !------------------------------------------------------------
  !NEW orimo
  complex(c_double_complex) :: charge
  integer(c_long) :: dimbc
  dimbc = min(dim, irad_ecs-1)
  !------------------------------------------------------------

  ld = ndvr + 1
  allocate(rrho2(1:dim, 1:2))

  do l = lll, ull
     !------------------------------------------------------------
     !NEW orimo
     bas_d2fac1(l) = one / xrad(dimbc + 1) ** (2 * l + 1)
     !------------------------------------------------------------
  end do


  do ifun = ifun_ll, ifun_ul
! Sato_ECS
     !------------------------------------------------------------
     !NEW orimo
     !do jfun = max(jfun_ll, ifun), max(jfun_ul, ifun)
     do jfun = jfun_ll, jfun_ul !! temporal code. this include some useless calculation
     !------------------------------------------------------------
!     do jfun = max(jfun_ll, ifun), max(jfun_ul, ifun)
! Sato_ECS
        do l = lll, ull
           m = -mval(ifun) + mval(jfun)
           if (l < abs(m)) then
              v2sph(1:dim, l, jfun, ifun) = czero
              cycle
           end if
           ! r.h.s and BC
           tmpr = zero
           tmpi = zero
           !------------------------------------------------------------
           !NEW orimo
           do irad = 1, dimbc
           !------------------------------------------------------------
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
           !------------------------------------------------------------
           !NEW orimo
           call DPBTRS('L', dimbc, ndvr, 2, d2ll(1,1,l), ld, rrho2, dim, info)
           !------------------------------------------------------------
           if (info /= 0) then
              write(6, "('hprod_poisson: dpbtrs bad info. ', 5i5)") ifun, jfun, l, m, info
              stop
           end if
           ! solution with the correct BC
           !------------------------------------------------------------
           !NEW orimo
           do irad = 1, dimbc
              v2sph(irad, l, jfun, ifun) = (rrho2(irad, 1) + tmpr * bas_d2rpl1(irad, l)) * bas_d2fac2(irad, l) &
                                         + (rrho2(irad, 2) + tmpi * bas_d2rpl1(irad, l)) * bas_d2fac2(irad, l) * iunit
           end do
           charge = (tmpr + tmpi * iunit) / bas_d2fac1(l)
           do irad = dimbc + 1, dim
              !! mean field in the scaled region
              v2sph(irad, l, jfun, ifun) = bas_d2crpl1(irad, l) * charge
           end do
        end do
     end do
  end do

  deallocate(rrho2)

end subroutine hprod_mkv2mf_poissonp_ecs
!######################################################################
subroutine hprod_mkv2mf_herm_ecs(v2_type, v2sph, v2ang)

  use, intrinsic :: iso_c_binding
  use mod_sph, only : lmax2, nlat
  use mod_rad, only : nrad, nradfc
  use mod_ormas, only : nfcore, nfcore1, nfcore2, nfun

  implicit none
  character(len=*), intent(in) :: v2_type
  complex(c_double_complex), intent(inout) :: v2sph(1:(nrad-1), 0:lmax2, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(inout) :: v2ang(1:(nrad-1), 1:nlat, 1:nfun, 1:nfun)

  integer(c_long) :: l, ilat, ifun, jfun, irad
  integer(c_long) :: dim, ifun_ll, ifun_ul, jfun_ll, jfun_ul, llr, ulr

  if (trim(v2_type) == 'tot') then
     dim = nrad - 1
     ifun_ll = 1
     ifun_ul = nfun
     jfun_ll = 1
     jfun_ul = nfun
  else if (trim(v2_type) == 'fcfc') then
!BUG dim = nradfc
     dim = nrad - 1
     ifun_ll = 1
     ifun_ul = nfcore
     jfun_ll = 1
     jfun_ul = nfcore
  else if (trim(v2_type) == 'fc1fc1') then
!BUG dim = nradfc
     dim = nrad - 1
     ifun_ll = nfcore2 + 1
     ifun_ul = nfcore
     jfun_ll = nfcore2 + 1
     jfun_ul = nfcore
  else if (trim(v2_type) == 'fc2fc2') then
!BUG dim = nradfc
     dim = nrad - 1
     ifun_ll = 1
     ifun_ul = nfcore2
     jfun_ll = 1
     jfun_ul = nfcore2
  else if (trim(v2_type) == 'fc1dyn') then
     dim = nradfc
     ifun_ll = nfcore2 + 1
     ifun_ul = nfcore
     jfun_ll = nfcore + 1
     jfun_ul = nfun
  else if (trim(v2_type) == 'fc2dyn') then
     dim = nradfc
     ifun_ll = 1
     ifun_ul = nfcore2
     jfun_ll = nfcore + 1
     jfun_ul = nfun
  else if (trim(v2_type) == 'dyn') then
     dim = nrad - 1
     ifun_ll = nfcore + 1
     ifun_ul = nfun
     jfun_ll = nfcore + 1
     jfun_ul = nfun
  else
     write(6, "('v2_type: ', 6A)") trim(v2_type)
     stop "bad v2_type in hprod_mkv2mf_herm"
  end if

  !$omp parallel default(shared) private(llr, ulr)
  !###########################
  call util_omp_disp(1, dim, llr, ulr)
  do ifun = ifun_ll, ifun_ul
     do jfun = max(jfun_ll, ifun + 1), jfun_ul
        do l = 0, lmax2
           do irad = llr, ulr
              v2sph(irad, l, ifun, jfun) = conjg(v2sph(irad, l, jfun, ifun))
           end do
        end do
        do ilat = 1, nlat
           do irad = llr, ulr
              v2ang(irad, ilat, ifun, jfun) = conjg(v2ang(irad, ilat, jfun, ifun))
           end do
        end do
     end do
  end do
  !###########################
  !$omp end parallel

end subroutine hprod_mkv2mf_herm_ecs
!######################################################################
