!#######################################################################
subroutine hprod_init()

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad
  use mod_sph, only : nlat, lmax1, lmax2
  use mod_control, only : istdcis, tdcis_rvg
  use mod_ormas, only : nact, nfun, lcic
  use mod_hprod, only : orb, torb, h0orb, h1orb, v2orb, gorb
  use mod_hprod, only : orbg, torbg, h0orbg, h1orbg, v2orbg, gorbg
  use mod_hprod, only : orb0, rhofc, v2jfc, v2xfc, rho1, rho2, v2sph, v2ang, cic, dcic, cic0
  use mod_hprod, only : ovlp, fmat, bmat, xmat, int1e, int2e, den1, rden, rrden, den2, den2r
  use mod_hprod, only : orbe, gorbe, v2orbe, v2ange, orbo, gorbo, v2orbo, v2ango, rho23j, v2sph3j, torb3j
! tdcis-teramura
  use mod_hprod, only : tdcis_eig, h1orb0, orb0rot, orb0rotg, v2ang0, &
       ridm_tdcis, ezphi, h0orb0
! tdcis-teramura

! Sato_ECS
  use mod_rad, only : ecs_flag
  use mod_hprod, only : h0orb_out, h1orb_out
! Sato_ECS

  implicit none
  integer(c_int) :: numx

  allocate(orb(1:(nrad-1),0:lmax1,1:nfun))
  allocate(torb(1:(nrad-1),0:lmax1,1:nfun))
  allocate(h0orb(1:(nrad-1),0:lmax1,1:nfun))
  allocate(h1orb(1:(nrad-1),0:lmax1,1:nfun))
  allocate(v2orb(1:(nrad-1),0:lmax1,1:nfun))
  allocate(gorb(1:(nrad-1),0:lmax1,1:nfun))
  allocate(orbg(1:(nrad-1),1:nlat,1:nfun))
  allocate(torbg(1:(nrad-1),1:nlat,1:nfun))
  allocate(h0orbg(1:(nrad-1),1:nlat,1:nfun))
  allocate(h1orbg(1:(nrad-1),1:nlat,1:nfun))
  allocate(v2orbg(1:(nrad-1),1:nlat,1:nfun))
  allocate(gorbg(1:(nrad-1),1:nlat,1:nfun))
  allocate(orb0(1:(nrad-1),0:lmax1,1:nfun))
  allocate(rhofc(1:(nrad-1)))
  allocate(v2jfc(1:(nrad-1)))
  allocate(v2xfc(1:(nrad-1)))
  allocate(rho1(1:(nrad-1),1:nlat))
  allocate(rho2(1:(nrad-1),0:lmax2,1:nfun,1:nfun))
  allocate(v2sph(1:(nrad-1),0:lmax2,1:nfun,1:nfun))
  allocate(v2ang(1:(nrad-1),1:nlat, 1:nfun,1:nfun))
  allocate(cic(1:lcic))
  allocate(dcic(1:lcic))
  allocate(cic0(1:lcic))

! tdcis-teramura
  if (istdcis) then
     allocate(tdcis_eig(1:nfun))
     allocate(h0orb0(1:(nrad-1),0:lmax1,1:nfun))
     allocate(h1orb0(1:(nrad-1),0:lmax1,1:nfun))
     allocate(orb0rot(1:(nrad-1),0:lmax1,1:nfun))
     allocate(orb0rotg(1:(nrad-1),1:nlat,1:nfun))
     allocate(v2ang0(1:(nrad-1),1:nlat, 1:nfun,1:nfun))
     allocate(ridm_tdcis(1:nfun,1:nfun))
     if (tdcis_rvg) then 
        allocate(ezphi(1:(nrad-1),0:lmax1,1:nfun))
     end if
  end if
! tdcis-teramura

! Sato_ECS
  if (ecs_flag == 1) then
     allocate(h0orb_out(1:(nrad-1),0:lmax1,1:nfun))
     allocate(h1orb_out(1:(nrad-1),0:lmax1,1:nfun))
  end if
! Sato_ECS

! ##### 3j selection rule #####
  !if (exact3j) then
  !   allocate(orbe  (1:(nrad-1),1:nlat,1:nfun))
  !   allocate(gorbe (1:(nrad-1),1:nlat,1:nfun))
  !   allocate(v2orbe(1:(nrad-1),1:nlat,1:nfun))
  !   allocate(v2ange(1:(nrad-1),1:nlat,1:nfun,1:nfun))
  !   allocate(orbo  (1:(nrad-1),1:nlat,1:nfun))
  !   allocate(gorbo (1:(nrad-1),1:nlat,1:nfun))
  !   allocate(v2orbo(1:(nrad-1),1:nlat,1:nfun))
  !   allocate(v2ango(1:(nrad-1),1:nlat,1:nfun,1:nfun))
  !   allocate(rho23j(1:(nrad-1),0:lmax2,1:nfun,1:nfun))
  !   allocate(v2sph3j(1:(nrad-1),0:lmax2,1:nfun,1:nfun))
  !   allocate(torb3j(1:(nrad-1),1:nlat,1:nfun,1:2))
  !end if
! ##### 3j selection rule #####

  allocate(ovlp(1:nfun,1:nfun))
  allocate(fmat(1:nfun,1:nfun))
  allocate(bmat(1:nfun,1:nfun))
  allocate(xmat(1:nfun,1:nfun))

  numx = max(1, nact)
  allocate(int1e(1:numx,1:numx))
  allocate(int2e(1:numx,1:numx,1:numx,1:numx))
  allocate(den1(1:numx,1:numx))
  allocate(rden(1:numx,1:numx))
  allocate(rrden(1:numx,1:numx))
  allocate(den2(1:numx,1:numx,1:numx,1:numx))
  allocate(den2r(1:numx,1:numx,1:numx,1:numx))

end subroutine hprod_init
!#######################################################################
