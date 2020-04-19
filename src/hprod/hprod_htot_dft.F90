!#######################################################################
subroutine hprod_htot_dft(dtime, lfield, wfn, hwfn)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : nbas
  use mod_ormas, only : nfun
  use mod_hprod, only : orb,h0orb,h1orb,gorb,v2orb

  implicit none
  real(c_double), intent(in) :: dtime, lfield(1:3,1:3)
  complex(c_double_complex), intent(in) :: wfn(1:nbas,1:nfun)
  complex(c_double_complex), intent(out) :: hwfn(1:nbas,1:nfun)

  call hprod_htotx_dft(dtime, lfield, wfn)
  hwfn(1:nbas,1:nfun) = 0d0
  call hprod_htot_dtorb(dtime, orb, h0orb, h1orb, gorb, v2orb, hwfn)

end subroutine hprod_htot_dft
!#######################################################################
subroutine hprod_htotx_dft(dtime, lfield, wfn)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : nbas, nbas2, ngrid
  use mod_control, only : ioorot, isplit, iprojfc, igauge, PSP
  use mod_const, only : zero, czero, runit, ctrue
  use mod_ormas, only : neltot, nfun, nfcore
  use mod_hprod, only : rho2, v2sph, v2ang
  use mod_hprod, only : orb, orbg, h0orb, h1orb, gorb, gorbg, v2orb
  use mod_hprod, only : ene_fcore, ene_dcore, ene_core, ene_act, ene_tot

  implicit none
  real(c_double), intent(in) :: dtime, lfield(1:3, 1:3)
  complex(c_double_complex), intent(in) :: wfn(1:nbas, 1:nfun)

  real(c_double) :: ene0, ene1, edftj, edftx
  real(c_double), external :: hprod_ene_dft1
  complex(c_double_complex) :: zfield, efield, afield
  complex(c_double_complex), allocatable :: rhoj(:)
  complex(c_double_complex), allocatable :: v2js(:)
  complex(c_double_complex), allocatable :: v2jg(:)
  real(c_double), allocatable :: rhoxc(:)

  allocate(rhoj(1:nbas2))
  allocate(v2js(1:nbas2))
  allocate(v2jg(1:ngrid))
  allocate(rhoxc(1:ngrid))

  ! laser field
  zfield = lfield(3, 1) ! chosen gauge
  efield = lfield(3, 2) ! electric field
  afield = lfield(3, 3) ! vector potential

  ! scratch initialization
  call hprod_htot_init

  ! current orbitals
  call hprod_orbin(lfield, wfn, orb, orbg)

  ! atomic hamiltonian
  if (ioorot == 0 .or. isplit == 0) then
     call hprod_tprod_dyn(orb, h0orb)
  end if
  if (iprojfc == 2) then
     call hprod_tprod_fc(orb, h0orb);
  end if
  if (PSP) call hprod_projpp(runit, orb, h0orb)

  ! laser hamiltonian
  if (ioorot .ne. 1 .or. isplit .ne. 1) then
    if (igauge == 0) then
       call hprod_zprod_dyn(zfield, orb, h1orb)
    else
       call hprod_pzprod_dyn(zfield, orb, h1orb)
    end if
  end if
  if (iprojfc == 2) then
     if (igauge == 0) then
        call hprod_zprod_fc(zfield, orb, h1orb)
     else
        call hprod_pzprod_fc(zfield, orb, h1orb)
     end if
  end if
  !NEW
  if (nfcore > 0 .and. igauge == 1) call hprod_zprod_fc(efield, orb, h1orb);
  !NEW

  ene0 = hprod_ene_dft1(orb, h0orb)
  ene1 = hprod_ene_dft1(orb, h1orb)

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

  ene_tot = ene_fcore + ene0 + ene1 + edftj + edftx
  call hprod_htot_mkxmat_aa1(dtime, orb, h0orb, h1orb, gorb, v2orb)

  deallocate(rhoxc)
  deallocate(v2jg)
  deallocate(v2js)
  deallocate(rhoj)

end subroutine hprod_htotx_dft
!#######################################################################
