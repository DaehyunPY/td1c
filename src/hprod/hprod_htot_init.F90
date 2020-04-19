!#######################################################################
subroutine hprod_htot_init()

  use, intrinsic :: iso_c_binding
  use mod_rad, only : ecs_flag
  use mod_bas, only : nbas, nbas2, ngrid
  use mod_ormas, only : nfun, nact, lcic
  use mod_hprod

  implicit none
  integer(c_long) :: size_orb, size_orbg
  integer(c_long) :: size_cic, size_den1, size_den2
  integer(c_long) :: size_fmat, size_rho2, size_rho2g

  ! initialization
  size_orb = nbas * nfun
  size_orbg = ngrid * nfun
  size_cic = lcic
  size_den1 = nact * nact
  size_den2 = size_den1 * size_den1
  size_fmat = nfun * nfun
  size_rho2 = nbas2 * nfun * nfun
  size_rho2g = ngrid * nfun * nfun

  call zclear_omp(size_cic, dcic)

  call zclear_omp(size_den1, den1)
  call zclear_omp(size_den1, rden)
  call zclear_omp(size_den1, rrden)
  call zclear_omp(size_den1, int1e)

  call zclear_omp(size_den2, den2)
  call zclear_omp(size_den2, int2e)

  call zclear_omp(size_fmat, fmat)
  call zclear_omp(size_fmat, bmat)
  call zclear_omp(size_fmat, xmat)

  call zclear_omp(size_orb, orb)  ; call zclear_omp(size_orbg, orbg)
  call zclear_omp(size_orb, torb) ; call zclear_omp(size_orbg, torbg)
  call zclear_omp(size_orb, h0orb); call zclear_omp(size_orbg, h0orbg)
  call zclear_omp(size_orb, h1orb); call zclear_omp(size_orbg, h1orbg)
  call zclear_omp(size_orb, gorb) ; call zclear_omp(size_orbg, gorbg)
  call zclear_omp(size_orb, v2orb); call zclear_omp(size_orbg, v2orbg)
  !##### 3j selection rule #####
  !oldif (exact3j) then
  !old   call zclear_omp(size_orbg, orbe);   call zclear_omp(size_orbg, orbo) 
  !old   call zclear_omp(size_orbg, gorbe);  call zclear_omp(size_orbg, gorbo)
  !old   call zclear_omp(size_orbg, v2orbe); call zclear_omp(size_orbg, v2orbo)
  !oldend if
  !##### 3j selection rule #####

! Sato_ECS
  if (ecs_flag == 1) then
     h0orb_out = 0d0
     h1orb_out = 0d0
  end if
! Sato_ECS

  ! carefull initialization for 2d quantities
  call hprod_htot_init2()

end subroutine hprod_htot_init
!#######################################################################
subroutine hprod_htot_init2()

  use, intrinsic :: iso_c_binding
  use mod_const, only : czero
  use mod_rad, only : nrad
  use mod_sph, only : nlat, lmax2
  use mod_ormas, only : nfcore, nfun
  use mod_hprod, only : rho2, v2sph, v2ang, v2ange, v2ango, rho23j, v2sph3j

  implicit none
  integer(c_long) :: ifun, jfun

  do ifun = 1, nfun
     do jfun = max(nfcore + 1, ifun), nfun
        rho2(1:(nrad-1), 0:lmax2, jfun, ifun) = czero
        rho2(1:(nrad-1), 0:lmax2, ifun, jfun) = czero
        v2sph(1:(nrad-1), 0:lmax2, jfun, ifun) = czero
        v2sph(1:(nrad-1), 0:lmax2, ifun, jfun) = czero
        v2ang(1:(nrad-1), 1:nlat, jfun, ifun) = czero
        v2ang(1:(nrad-1), 1:nlat, ifun, jfun) = czero
        !##### 3j selection rule #####
        !oldif (exact3j) then
        !old   v2ange(1:(nrad-1), 1:nlat, jfun, ifun) = czero
        !old   v2ange(1:(nrad-1), 1:nlat, ifun, jfun) = czero
        !old   v2ango(1:(nrad-1), 1:nlat, jfun, ifun) = czero
        !old   v2ango(1:(nrad-1), 1:nlat, ifun, jfun) = czero
        !old   rho23j(1:(nrad-1), 0:lmax2, jfun, ifun) = czero
        !old   rho23j(1:(nrad-1), 0:lmax2, ifun, jfun) = czero
        !old   v2sph3j(1:(nrad-1), 0:lmax2, jfun, ifun) = czero
        !old   v2sph3j(1:(nrad-1), 0:lmax2, ifun, jfun) = czero
        !oldend if
        !##### 3j selection rule #####
     end do
  end do

end subroutine hprod_htot_init2
!#######################################################################
