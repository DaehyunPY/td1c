!#######################################################################
subroutine hprod_htot(dtime, lfield, wfn, cic, hwfn, hcic)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : nbas
  use mod_control, only : cionly
  use mod_ormas, only : lcic,nfun
  use mod_hprod, only : dcic,orb,h0orb,h1orb,gorb,v2orb
  use mod_control, only : istdcis

  implicit none
  real(c_double), intent(in) :: dtime, lfield(1:3, 1:3)
  complex(c_double_complex), intent(in) :: wfn(1:*)
  complex(c_double_complex), intent(in) :: cic(1:lcic)
  complex(c_double_complex), intent(out) :: hwfn(1:*)
  complex(c_double_complex), intent(out) :: hcic(1:lcic)

! tdcis-teramura
  if (istdcis) then
     call hprod_htot_tdcis(dtime, lfield, wfn, cic, hwfn, hcic)
     return
  end if
! tdcis-teramura

  call hprod_htotx(dtime, lfield, wfn, cic)
  hcic(1:lcic) = dcic(1:lcic)
  hwfn(1:nbas*nfun) = 0d0
  if (.not.cionly) call hprod_htot_dtorb(dtime, orb, h0orb, h1orb, gorb, v2orb, hwfn)

!  if (projhigh) then
!     call hprod_projhigh(hwfn)
!  end if

!debug
!  write(6, "('hprod_htot: wfn')")
!  call hprod_printorb(wfn)
!  write(6, "('hprod_htot: hwfn')")
!  call hprod_printorb(hwfn)
!  stop '@ hprod_htot for debug.'
!debug

end subroutine hprod_htot
!#######################################################################
subroutine hprod_htotx(dtime, lfield, wfn, cic)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : ecs_flag
  use mod_control, only : exact3j
  use mod_bas, only : nbas, ngrid
  use mod_const, only : zero, czero, runit, cfalse
  use mod_hprod, only : rho1, rho2, v2sph, v2ang, v2ange, v2ango
  use mod_ormas, only : neltot, nfun, nfcore, ndcore, nact, lcic
  use mod_control, only : ioorot, isplit, iprojfc, igauge, dft_type, PSP, istdcis
  use mod_hprod, only : ene_fcore, ene_dcore, ene_core, ene_act, ene_tot
  use mod_hprod, only : dcic, den1, den2, rden, rrden, den2r, int1e, int2e
  use mod_hprod, only : orb, orbg, h0orb, h1orb, gorb, gorbg, v2orb, v2orbg
  use mod_hprod, only : gorbe,gorbo,v2orbe,v2orbo,torb3j,rho23j,v2sph3j,orbe,orbo

  implicit none
  real(c_double), intent(in) :: dtime, lfield(1:3, 1:3)
  complex(c_double_complex), intent(in) :: wfn(1:*)
  complex(c_double_complex), intent(in) :: cic(1:lcic)

  integer(c_int) :: size_orb, size_forb, size_cic
  real(c_double), external :: hprod_ene_fcore
  real(c_double), external :: hprod_ene_dcore
  real(c_double), external :: hprod_ene_act
  complex(c_double_complex) :: zfield, efield, afield

  if (dft_type .ne. 0) then
     dcic(1:lcic) = czero
     call hprod_htotx_dft(dtime, lfield, wfn)
     return
! tdcis-teramura
  else if (istdcis) then !!teramura
     call hprod_htotx_tdcis(dtime, lfield, wfn, cic)
     return
! tdcis-teramura
  end if

  ! laser field
  zfield = lfield(3, 1) ! chosen gauge
  efield = lfield(3, 2) ! electric field
  afield = lfield(3, 3) ! vector potential

  ! scratch initialization
  call hprod_htot_init
  ! current orbitals
  call hprod_orbin(lfield, wfn, orb, orbg)

  ! RDMs
  call ormas_mkden1(cic, den1)
  call ormas_mkden2(cic, den1, den2)
  call hprod_invden(den1, rden, rrden)

  ! atomic hamiltonian
  if (ioorot == 0 .or. isplit == 0) then
     call hprod_tprod_dyn(orb, h0orb)
  end if
  if (iprojfc == 2) then
     call hprod_tprod_fc(orb, h0orb);
  end if

  ! pseudopotential
  if (PSP) call hprod_projpp(runit, lfield, orb, h0orb)

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

  ! meanfield operator
  if (neltot(3) >= 2) then
     if (.not.exact3j) then
! Orimo_ECS
        if (ecs_flag == 0) then
           call hprod_mkrho2_dyn
           call hprod_mkv2mf_dyn
           !call hprod_mkv2mf2_dyn
        else
           call hprod_mkrho2_dyn_ecs
           call hprod_mkv2mf_dyn_ecs
        end if
! Orimo_ECS
        !call hprod_mfprod
        call hprod_mfprod2
        call bas_ang2sph1_dyn(gorbg, gorb);
        call bas_ang2sph1_dyn(v2orbg, v2orb);
     else
        call hprod_mkrho2_x3j_dyn
        call hprod_mkv2mf_x3j_dyn
        call hprod_mfprodx3j
!old        ! ##### 3j selection rule #####
!old        call hprod_mfprod3j()
!old        call bas_ang2sph1_dyn3j(gorbe, gorbo, gorb);
!old        call bas_ang2sph1_dyn3j(v2orbe, v2orbo, v2orb);
!old        ! ##### 3j selection rule #####
     end if
  end if

!  if (projhigh) then
!     call hprod_projhigh(h0orb)
!     call hprod_projhigh(h1orb)
!     call hprod_projhigh(gorb)
!  end if

  ! mo integrals
  call hprod_mkint1_sph(cfalse, orb, h0orb, h1orb, gorb)
  if (.not.exact3j) then
     call hprod_mkint2_sph(rho2, v2sph)
  else
     call hprod_mkint2_x3j(rho2, v2sph)
!old     ! ##### 3j selection rule #####
!old     call hprod_mkint2_sph(rho23j, v2sph3j)
!old     ! ##### 3j selection rule #####
  end if

  ! energy
  ene_dcore = zero
  ene_core = zero
  ene_act = zero
  if (ndcore > 0) ene_dcore = hprod_ene_dcore(cfalse, orb, h0orb, h1orb, gorb)
  if (nact > 0) ene_act = hprod_ene_act(int1e, int2e, den1, den2)
  ene_core = ene_fcore + ene_dcore
  ene_tot = ene_core + ene_act
!DEBUG
!  write(6, "('hprod_htot: enefc  = ', f20.10)") ene_fcore
!  write(6, "('hprod_htot: enedc  = ', f20.10)") ene_dcore
!  write(6, "('hprod_htot: eneact = ', f20.10)") ene_act
!  write(6, "('hprod_htot: enetot = ', f20.10)") ene_tot
!  stop
!DEBUG

  ! sigma vector
  ! included in ormas_xact2ras
  ! call ormas_hcic(int1e, int2e, cic, dcic, ene_act)

! Orimo_ECS
  if (ecs_flag == 0) then
     call hprod_htot_mkxmat_fc(dtime, orb, h0orb, h1orb, gorb, v2orb)
     call hprod_htot_mkxmat_cc(dtime, orb, h0orb, h1orb, gorb, v2orb)
     call hprod_htot_mkxmat_ca(dtime, orb, h0orb, h1orb, gorb, v2orb)
     call hprod_htot_mkxmat_aa1(dtime, orb, h0orb, h1orb, gorb, v2orb)
     call hprod_htot_mkxmat_aa2(dtime, orb, h0orb, h1orb, gorb, v2orb, cic, dcic)
  else
     call hprod_htot_mkxmat_fc_ecs(dtime, orb, h0orb, h1orb, gorb, v2orb)
     call hprod_htot_mkxmat_cc_ecs(dtime, orb, h0orb, h1orb, gorb, v2orb)
     call hprod_htot_mkxmat_ca_ecs(dtime, orb, h0orb, h1orb, gorb, v2orb)
     call hprod_htot_mkxmat_aa1_ecs(dtime, orb, h0orb, h1orb, gorb, v2orb)
     call hprod_htot_mkxmat_aa2_ecs(dtime, orb, h0orb, h1orb, gorb, v2orb, cic, dcic)
  end if
! Orimo_ECS

end subroutine hprod_htotx
!#######################################################################
subroutine hprod_htoto(dtime, lfield, wfn, cic, hwfn)

  use, intrinsic :: iso_c_binding
  use mod_control, only : exact3j
  use mod_bas, only : nbas, ngrid
  use mod_const, only : zero, czero, runit, cfalse
  use mod_hprod, only : rho1, rho2, v2sph, v2ang, v2ange, v2ango
  use mod_ormas, only : neltot, nfun, nfcore, ndcore, nact, lcic
  use mod_control, only : ioorot, isplit, iprojfc, igauge, dft_type, PSP
  use mod_hprod, only : ene_fcore, ene_dcore, ene_core, ene_act, ene_tot
  use mod_hprod, only : den1, den2, rden, rrden, den2r, int1e, int2e
  use mod_hprod, only : orb, orbg, h0orb, h1orb, gorb, gorbg, v2orb, v2orbg
  use mod_hprod, only : gorbe,gorbo,v2orbe,v2orbo,torb3j,rho23j,v2sph3j,orbe,orbo

  implicit none
  real(c_double), intent(in) :: dtime, lfield(1:3, 1:3)
  complex(c_double_complex), intent(in) :: wfn(1:*)
  complex(c_double_complex), intent(in) :: cic(1:lcic)
  complex(c_double_complex), intent(out) :: hwfn(1:*)

  integer(c_int) :: size_orb, size_forb, size_cic
  real(c_double), external :: hprod_ene_fcore
  real(c_double), external :: hprod_ene_dcore
  real(c_double), external :: hprod_ene_act
  complex(c_double_complex) :: zfield, efield, afield

  stop 'hprod_htoto is nyi.'
  hwfn(1:nbas*nfun) = czero
  
  if (dft_type .ne. 0) then
     call hprod_htotx_dft(dtime, lfield, wfn)
     return
  end if
  ! laser field
  zfield = lfield(3, 1) ! chosen gauge
  efield = lfield(3, 2) ! electric field
  afield = lfield(3, 3) ! vector potential

  ! scratch initialization
  call hprod_htot_init

  ! current orbitals
  call hprod_orbin(lfield, wfn, orb, orbg)

  ! RDMs
  call ormas_mkden1(cic, den1)
  call ormas_mkden2(cic, den1, den2)
  call hprod_invden(den1, rden, rrden)

  ! atomic hamiltonian
  if (ioorot == 0 .or. isplit == 0) then
     call hprod_tprod_dyn(orb, h0orb)
  end if
  if (iprojfc == 2) then
     call hprod_tprod_fc(orb, h0orb);
  end if

  ! pseudopotential
  if (PSP) call hprod_projpp(runit, lfield, orb, h0orb)

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

  ! meanfield operator
  if (neltot(3) >= 2) then
     if (.not.exact3j) then
        call hprod_mkrho2_dyn
        call hprod_mkv2mf_dyn
        !call hprod_mkv2mf2_dyn
        !call hprod_mfprod
        call hprod_mfprod2
        call bas_ang2sph1_dyn(gorbg, gorb);
        call bas_ang2sph1_dyn(v2orbg, v2orb);
     else
        call hprod_mkrho2_x3j_dyn
        call hprod_mkv2mf_x3j_dyn
        call hprod_mfprodx3j
     end if
  end if

  ! mo integrals
  call hprod_mkint1_sph(cfalse, orb, h0orb, h1orb, gorb)
  if (.not.exact3j) then
     call hprod_mkint2_sph(rho2, v2sph)
  else
     call hprod_mkint2_x3j(rho2, v2sph)
!old     ! ##### 3j selection rule #####
!old     call hprod_mkint2_sph(rho23j, v2sph3j)
!old     ! ##### 3j selection rule #####
  end if

  ! energy
  ene_dcore = zero
  ene_core = zero
  ene_act = zero
  if (ndcore > 0) ene_dcore = hprod_ene_dcore(cfalse, orb, h0orb, h1orb, gorb)
  if (nact > 0) ene_act = hprod_ene_act(int1e, int2e, den1, den2)
  ene_core = ene_fcore + ene_dcore
  ene_tot = ene_core + ene_act

  ! orbital and ci derivatives
  size_orb = nbas * nfun
  size_forb = nbas * nfcore
  call hprod_htoto_getdt(dtime)
!  call zcopy_omp(size_orb, v2orb, hwfn)
  call zcopy(size_orb, v2orb, 1, hwfn, 1)
  call zclear_omp(size_forb, hwfn)

end subroutine hprod_htoto
!#######################################################################
