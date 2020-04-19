!#######################################################################
subroutine hprod_final()

  use, intrinsic :: iso_c_binding
  use mod_control, only : istdcis, tdcis_rvg
  use mod_hprod, only : orb, torb, h0orb, h1orb, v2orb, gorb
  use mod_hprod, only : orbg, torbg, h0orbg, h1orbg, v2orbg, gorbg
  use mod_hprod, only : orb0, rhofc, v2jfc, v2xfc, rho1, rho2, v2sph, v2ang, cic, dcic, cic0
  use mod_hprod, only : ovlp, fmat, bmat, xmat, int1e, int2e, den1, rden, rrden, den2, den2r
  use mod_hprod, only : orbe, gorbe, v2orbe, v2ange, orbo, gorbo, v2orbo, v2ango, rho23j, v2sph3j, torb3j
  use mod_hprod, only : projhigh, projhigh_orbs, projhigh_eigs, projhigh_ncut
! Sato_ECS
  use mod_rad, only : ecs_flag
  use mod_hprod, only : h0orb_out, h1orb_out
! Sato_ECS
  !!tdcis
  use mod_hprod, only : tdcis_eig, h1orb0, orb0rot, orb0rotg, &
       v2ang0, ridm_tdcis, ezphi, h0orb0
  !!tdcis

  implicit none

  deallocate(orb)
  deallocate(torb)
  deallocate(h0orb)
  deallocate(h1orb)
  deallocate(v2orb)
  deallocate(gorb)
  deallocate(orbg)
  deallocate(torbg)
  deallocate(h0orbg)
  deallocate(h1orbg)
  deallocate(v2orbg)
  deallocate(gorbg)
  deallocate(orb0)
  deallocate(rhofc)
  deallocate(v2jfc)
  deallocate(v2xfc)
  deallocate(rho1)
  deallocate(rho2)
  deallocate(v2sph)
  deallocate(v2ang)
  deallocate(cic)
  deallocate(dcic)
  deallocate(cic0)

! Sato_ECS
  if (ecs_flag == 1) then
     deallocate(h0orb_out)
     deallocate(h1orb_out)
  end if
! Sato_ECS

! tdcis-teramura
  if(istdcis) then
     deallocate(tdcis_eig)
     deallocate(h0orb0)
     deallocate(h1orb0)
     deallocate(orb0rot)
     deallocate(orb0rotg)
     deallocate(v2ang0)
     deallocate(ridm_tdcis)
     if(tdcis_rvg) then
        deallocate(ezphi)
     end if
  end if
! tdcis-teramura

  if (projhigh) then
     deallocate(projhigh_orbs)
     deallocate(projhigh_eigs)
     deallocate(projhigh_ncut)
  end if

! ##### 3j selection rule #####
  !if (exact3j) then
  !   deallocate(orbe)
  !   deallocate(gorbe)
  !   deallocate(v2orbe)
  !   deallocate(v2ange)
  !   deallocate(orbo)
  !   deallocate(gorbo)
  !   deallocate(v2orbo)
  !   deallocate(v2ango)
  !   deallocate(rho23j)
  !   deallocate(v2sph3j)
  !   deallocate(torb3j)
  !end if
! ##### 3j selection rule #####

  deallocate(ovlp)
  deallocate(fmat)
  deallocate(bmat)
  deallocate(xmat)

  deallocate(int1e)
  deallocate(int2e)
  deallocate(den1)
  deallocate(rden)
  deallocate(rrden)
  deallocate(den2)
  deallocate(den2r)

end subroutine hprod_final
!#######################################################################
