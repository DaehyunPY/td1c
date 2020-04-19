!#######################################################################
subroutine hprod_energy_dft(lfield, wfn)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : nbas2, ngrid
  use mod_control, only : igauge, PSP
  use mod_const, only : zero, czero, runit, ctrue
  use mod_ormas, only : neltot, nfun, nfcore
  use mod_hprod, only : rho2, v2sph, v2ang
  use mod_hprod, only : orb, orbg, h1orb, gorb, gorbg
  use mod_hprod, only : ene_fcore, ene_tot

  implicit none
  real(c_double), intent(in) :: lfield(1:3, 1:3)
  complex(c_double_complex), intent(in) :: wfn(*)

  real(c_double) :: ene1, edftj, edftx
  real(c_double), external :: hprod_ene_dft1
  complex(c_double_complex) :: zfield
  complex(c_double_complex), allocatable :: rhoj(:)
  complex(c_double_complex), allocatable :: v2js(:)
  complex(c_double_complex), allocatable :: v2jg(:)
  real(c_double), allocatable :: rhoxc(:)

  allocate(rhoj(1:nbas2))
  allocate(v2js(1:nbas2))
  allocate(v2jg(1:ngrid))
  allocate(rhoxc(1:ngrid))

  ! laser field
  zfield = lfield(3, 1)

  ! scratch initialization
  call hprod_htot_init

  ! current orbitals
  call hprod_orbin(lfield, wfn, orb, orbg)

  ! 1e operators
  call hprod_tprod_dyn(orb, h1orb)
  if (PSP) call hprod_projpp(runit, lfield, orb, h1orb)

  if (igauge == 0) call hprod_zprod_dyn(zfield, orb, h1orb)
  if (igauge == 1) call hprod_pzprod_dyn(zfield, orb, h1orb)
  ene1 = ene_fcore + hprod_ene_dft1(orb, h1orb)

!  call hprod_mkrho1(rhoxc)
  call hprod_mkrhoxc(rhoj, rhoxc)
  call hprod_poisson1(rhoj, v2js)
  call bas_sph2ang2one(0, v2js, v2jg)
  call hprod_mfprod_dftj(v2jg, orbg, gorbg, edftj)
  call hprod_mfprod_dftx(rhoxc, orbg, gorbg, edftx)
!scf  call bas_ang2sph1_dyn(gorbg, gorb);

  ene_tot = ene1 + edftj + edftx
  !debug
  write(6, "('hprod_energy_dft:', 4f20.10)") ene1, edftj, edftx, ene_tot
  !debug

!nyi  ene_dcore = zero
!nyi  ene_core = zero
!nyi  ene_act = zero
!nyi  if (ndcore > 0) ene_dcore = hprod_ene_dcore(ctrue, orb, h0orb, h1orb, gorb)
!nyi  if (nact > 0) ene_act = hprod_ene_act(int1e, int2e, den1, den2)
!nyi  ene_core = ene_fcore + ene_dcore
!nyi  ene_tot = ene_core + ene_act

  deallocate(rhoxc)
  deallocate(v2jg)
  deallocate(v2js)
  deallocate(rhoj)

end subroutine hprod_energy_dft
!#######################################################################
