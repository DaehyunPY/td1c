!#######################################################################
subroutine hprod_mpfock_dft(lfield, wfn)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : nbas, nbas2, ngrid
  use mod_control, only : ioorot, isplit, iprojfc, igauge, PSP
  use mod_const, only : zero, czero, runit, ctrue
  use mod_ormas, only : neltot, nfun, nfcore
  use mod_hprod, only : rho2, v2sph, v2ang
  use mod_hprod, only : orb, orbg, h0orb, h1orb, gorb, gorbg, v2orb, v2orbg
  use mod_hprod, only : ene_fcore, ene_dcore, ene_core, ene_act, ene_tot

  implicit none
  real(c_double), intent(in) :: lfield(1:3, 1:3)
  complex(c_double_complex), intent(in) :: wfn(1:nbas, 1:nfun)

  real(c_double) :: ene0, ene1, edftj, edftx
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

  call hprod_htot_init
  call hprod_orbin(lfield, wfn, orb, orbg)

  ! 1e operators
  zfield = lfield(3, 1)
  call hprod_tprod_all(orb, h0orb)
  if (PSP) call hprod_projpp(runit, lfield, orb, h0orb)

  if (igauge == 0) call hprod_zprod_all(zfield, orb, h1orb)
  if (igauge == 1) call hprod_pzprod_all(zfield, orb, h1orb)

!  call hprod_mkrho1(rhoxc)
  call hprod_mkrhoxc(rhoj, rhoxc)
  call hprod_poisson1(rhoj, v2js)
  call bas_sph2ang2one(0, v2js, v2jg)
  call hprod_mfprod_dftj(v2jg, orbg, gorbg, edftj)
  call hprod_mfprod_dftx(rhoxc, orbg, gorbg, edftx)
  call bas_ang2sph1_dyn(gorbg, gorb);
  if (nfcore > 0) then
     call bas_ang2sph1_fc(gorbg, gorb);
  end if

  !debug
  ene0 = hprod_ene_dft1(orb, h0orb)
  ene1 = hprod_ene_dft1(orb, h1orb)
  ene_tot = ene0 + ene1 + edftj + edftx
  write(6, "('hprod_mpfock_dft:', 5f20.10)") ene0, ene1, edftj, edftx, ene_tot
  !debug

  ! gfock vector
  call hprod_mpfock_sum()

!debug
!  call zscal(ngrid*nfun, czero, gorbg, 1)
!  call zscal(ngrid*nfun, czero, v2orbg, 1)
!  call hprod_mfprod_dftj(v2jg, orbg, gorbg, edftj)
!  call hprod_mfprod_dftx(rhoxc, orbg, v2orbg, edftx)
!  call bas_ang2sph1_dyn(gorbg, gorb);
!  call bas_ang2sph1_dyn(v2orbg, v2orb);
!  write(6, "('hprod_mpfock_dft:', 5f20.10)") ene0, ene1, edftj, edftx, ene_tot
!  edftj = hprod_ene_dft1(orb, gorb) * 0.5d+0
!  edftx = hprod_ene_dft1(orb, v2orb) * 3.d+0 / 4.d+0
!  ene_tot = ene0 + ene1 + edftj + edftx
!  write(6, "('hprod_mpfock_dft:', 5f20.10)") ene0, ene1, edftj, edftx, ene_tot
!  ene_tot = (ene0 + ene1) * 0.5d+0 + edftj + edftx * 4.d+0 / 3.d+0 / 2.d+0
!  write(6, "('hprod_mpfock_dft:', 5f20.10)") ene0, ene1, edftj, edftx, ene_tot
!  call zaxpy(nbas*nfun, runit, gorb, 1, h0orb, 1)
!  call zaxpy(nbas*nfun, runit, v2orb, 1, h0orb, 1)
!  ene_tot = hprod_ene_dft1(orb, h0orb) * 0.5d+0
!  write(6, "('hprod_mpfock_dft:', 5f20.10)") ene0, ene1, edftj, edftx, ene_tot
!  ene_tot = hprod_ene_dft1(orb, gorb) * 0.5d+0
!  write(6, "('hprod_mpfock_dft:', 5f20.10)") ene0, ene1, edftj, edftx, ene_tot
!debug

  deallocate(rhoxc)
  deallocate(v2jg)
  deallocate(v2js)
  deallocate(rhoj)

end subroutine hprod_mpfock_dft
!#######################################################################
