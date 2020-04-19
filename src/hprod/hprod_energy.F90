!#######################################################################
subroutine hprod_energy(lfield, wfn, cic)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : nbas, ngrid
  use mod_control, only : igauge, dft_type, PSP, exact3j
  use mod_const, only : zero, czero, runit, ctrue
  use mod_ormas, only : neltot, nfun, nfcore, ndcore, nact
  use mod_hprod, only : rho2, v2sph, v2ang
  use mod_hprod, only : orb, orbg, h0orb, h1orb, gorb, gorbg
  use mod_hprod, only : ene_fcore, ene_dcore, ene_core, ene_act, ene_tot
  use mod_hprod, only : dcic, den1, den2, rden, rrden, den2r, int1e, int2e
! Orimo_ECS
  use mod_rad, only : ecs_flag
! Orimo_ECS

  implicit none
  real(c_double), intent(in) :: lfield(1:3, 1:3)
  complex(c_double_complex), intent(in) :: wfn(*)
  complex(c_double_complex), intent(in) :: cic(*)

  integer(c_int) :: size_orb, size_forb, size_cic
  real(c_double), external :: hprod_ene_dcore
  real(c_double), external :: hprod_ene_act
  complex(c_double_complex) :: zfield
!debug
  integer(c_int) :: iact,jact
!debug

  if (dft_type .ne. 0) then
     call hprod_energy_dft(lfield, wfn)
     return
  end if

  ! laser field
  zfield = lfield(3, 1)

  ! scratch initialization
  call hprod_htot_init

  ! current orbitals
  call hprod_orbin(lfield, wfn, orb, orbg)

  ! RDMs
  call ormas_mkden1(cic, den1)
  call ormas_mkden2(cic, den1, den2)

  ! 1e operators
!New_ene_FC
!  call hprod_tprod_ene(orb, h0orb)
!  if (igauge == 0) call hprod_zprod_all(zfield, orb, h1orb)
!  if (igauge == 1) call hprod_pzprod_all(zfield, orb, h1orb)
  call hprod_tprod_dyn(orb, h0orb)
  if (PSP) call hprod_projpp(runit, lfield, orb, h0orb)

  if (igauge == 0) call hprod_zprod_dyn(zfield, orb, h1orb)
  if (igauge == 1) call hprod_pzprod_dyn(zfield, orb, h1orb)
!New_ene_FC

  ! meanfield operator
  if (neltot(3) >= 2) then
     if (.not.exact3j) then
! Orimo_ECS
        if (ecs_flag == 0) then
           call hprod_mkrho2_dyn
           call hprod_mkv2mf_dyn
        else
           call hprod_mkrho2_dyn_ecs
           call hprod_mkv2mf_dyn_ecs
        end if
! Orimo_ECS
        !call hprod_mkv2mf2_dyn
        !call hprod_mfprod_ene
        call hprod_mfprod2_ene
        call bas_ang2sph1_dyn(gorbg, gorb);
        !New_ene_FC
        !if (nfcore > 0) call bas_ang2sph1_fc(gorbg, gorb)
        !New_ene_FC
     else
        call hprod_mkrho2_x3j_dyn
        call hprod_mkv2mf_x3j_dyn
        call hprod_mfprodx3j_ene
     end if
  end if

  ! mo integrals
  call hprod_mkint1_sph(ctrue, orb, h0orb, h1orb, gorb)
  if (.not.exact3j) then
     call hprod_mkint2_sph(rho2, v2sph)
  else
     call hprod_mkint2_x3j(rho2, v2sph)
  end if
!DEBUG
!  call util_print_vec(nfun*nbas, orb, "test.orb")
!  call util_print_vec(nfun*nbas, h0orb, "test.h0orb")
!  call util_print_vec(nfun*nbas, h1orb, "test.h1orb")
!  call util_print_vec(nfun*nbas, gorb, "test.gorb")
!  call util_print_vec(nact**2, int1e, "test.int1e")
!  call util_print_vec(nact**4, int2e, "test.int2e")
!  stop "STOP for debug @ hprod_energy."
!DEBUG

  ! energy
  ene_dcore = zero
  ene_core = zero
  ene_act = zero
  if (ndcore > 0) ene_dcore = hprod_ene_dcore(ctrue, orb, h0orb, h1orb, gorb)
  if (nact > 0) ene_act = hprod_ene_act(int1e, int2e, den1, den2)
  ene_core = ene_fcore + ene_dcore
  ene_tot = ene_core + ene_act
!DEBUG
!  do iact = 1, nact
!     write(6,"('hprod_energy. int1e: ',i5,2e20.10)") iact,int1e(iact,iact)
!     do jact = 1, nact
!        write(6,"('hprod_energy. int2e: ',2i5,10e20.10)") iact,jact,int2e(iact,iact,jact,jact)
!     end do
!  end do
!  stop
!DEBUG

!DEBUG
!  write(6, "('hprod_energy: ene = ', 5f20.10)") ene_fcore, ene_dcore, ene_core, ene_act, ene_tot
!  stop
!DEBUG

end subroutine hprod_energy
!#######################################################################
real(c_double) function hprod_energy1(lfield, rion, wfn, cic)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : nbas, ngrid
  use mod_control, only : igauge, dft_type, PSP, exact3j
  use mod_const, only : zero, czero, runit, ctrue
  use mod_ormas, only : neltot, nfun, nfcore, ndcore, nact
  use mod_hprod, only : rho2, v2sph, v2ang
  use mod_hprod, only : orb, orbg, h0orb, h1orb, gorb, gorbg
  use mod_hprod, only : ene_fcore, ene_dcore, ene_core, ene_act, ene_tot
  use mod_hprod, only : int1e, int2e, den1, den2, ovlp

  implicit none
  real(c_double), intent(in) :: lfield(1:3, 1:3), rion
  complex(c_double_complex), intent(in) :: wfn(*)
  complex(c_double_complex), intent(in) :: cic(*)

  integer(c_int) :: size_orb, size_forb
  real(c_double), external :: hprod_ene_dcore
  real(c_double), external :: hprod_ene_act
  complex(c_double_complex) :: zfield

  if (dft_type .ne. 0) then
     call hprod_energy_dft(lfield, wfn)
     return
  end if

  ! laser field
  zfield = lfield(3, 1)

  !scratch initialization
  call hprod_htot_init

  ! current orbitals
  call hprod_orbin(lfield, wfn, orb, orbg)

  !!!!! special density !!!!!
  den1 = ovlp
  den2 = 0d0

  ! 1e operators
!New_ene_FC
!  call hprod_tprod_ene(orb, h0orb)
!  if (igauge == 0) call hprod_zprod_all(zfield, orb, h1orb)
!  if (igauge == 1) call hprod_pzprod_all(zfield, orb, h1orb)
  call hprod_tprod_dyn(orb, h0orb)
  if (PSP) call hprod_projpp(runit, lfield, orb, h0orb)

  if (igauge == 0) call hprod_zprod_dyn(zfield, orb, h1orb)
  if (igauge == 1) call hprod_pzprod_dyn(zfield, orb, h1orb)
!New_ene_FC

  ! mo integrals
  !call hprod_mkint1_sph(ctrue, orb, h0orb, h1orb, gorb)
  call hprod_mkint1_sph_inner(rion, orb, h0orb, h1orb, gorb)
  int2e = 0d0

  ! energy
  ene_dcore = zero
  ene_core = zero
  ene_act = zero
  if (ndcore > 0) ene_dcore = hprod_ene_dcore(ctrue, orb, h0orb, h1orb, gorb)
  if (nact > 0) ene_act = hprod_ene_act(int1e, int2e, den1, den2)
  ene_core = ene_fcore + ene_dcore
  ene_tot = ene_core + ene_act
  hprod_energy1 = ene_tot

!DEBUG
!  write(6, "('hprod_energy: ene = ', 5f20.10)") ene_fcore, ene_dcore, ene_core, ene_act, ene_tot
!  stop
!DEBUG

end function hprod_energy1
!#######################################################################
