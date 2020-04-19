!######################################################################
! test FC use nfcore_tdcis instead of nfcore
!#######################################################################
subroutine hprod_htot_tdcis(dtime, lfield, wfn, cic, hwfn, hcic)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : nbas
  use mod_ormas, only : lcic,nfun
  use mod_hprod, only : dcic,orb,h0orb,h1orb,gorb,v2orb
  use mod_hprod, only : projhigh
  use mod_control, only : tdcis_rvg

  !debug
  use mod_hprod, only : ovlp
  use mod_rad, only : xrad, nrad
  use mod_const, only : iunit
  use mod_sph, only : lmax1
  !debug
  
  implicit none
  real(c_double), intent(in) :: dtime
  real(c_double), intent(in) :: lfield(1:3, 1:3)
  complex(c_double_complex), intent(in) :: wfn(1:*)
  complex(c_double_complex), intent(in) :: cic(1:lcic)
  complex(c_double_complex), intent(out) :: hwfn(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(out) :: hcic(1:lcic)
  complex(c_double_complex) :: afield
  !debug
  ! integer(c_int) :: irad
  integer(c_int) :: ifun
  complex(c_double_complex) :: o_o(1:nfun,1:nfun)
  ! complex(c_double_complex) :: o_w(1:nfun,1:nfun)
  ! complex(c_double_complex) :: o_dto(1:nfun,1:nfun)
  ! complex(c_double_complex) :: o_h0(1:nfun,1:nfun)
  ! complex(c_double_complex) :: o_h1(1:nfun,1:nfun)
  ! complex(c_double_complex) :: dotorb(1:(nrad-1), 0:lmax1, 1:nfun)
  ! complex(c_double_complex) :: dotcic(1:nfun, 1:nfun)
  ! complex(c_double_complex) :: dotcic
  !debug

  afield = lfield(3, 3) ! vector potential

  ! call hprod_mkovlp(xrad(nrad), orb, orb, o_o)
  ! write(6 , "('hprod_htot_tdcis: debug_chichi ')") 
  ! do ifun = 1, nfun
  !    write(6 , "(2es20.10)") o_o(ifun, ifun)
  ! end do
  
  ! write(6 , "('hprod_htot_tdcis: debug_1 ')") 
  call hprod_htotx_tdcis(dtime, lfield, wfn, cic)
  hcic(1:lcic) = dcic(1:lcic)
  hwfn = 0d0

  ! write(6 , "('hprod_htot_tdcis: debug_2 ')") 
  if (tdcis_rvg) then
     call hprod_htot_dtorb_tdcis_rvg(dtime, orb, h0orb, h1orb, gorb, v2orb, cic, afield, hwfn)
  else 
     call hprod_htot_dtorb_tdcis(dtime, orb, h0orb, h1orb, gorb, v2orb, cic, hwfn)
  end if
  if (projhigh) then
     call hprod_projhigh(hwfn)
  end if

  !debug
  ! call hprod_mkovlp(xrad(nrad), orb, orb, o_o)
  ! call hprod_mkovlp(xrad(nrad), orb, hwfn, o_w)
  ! dotcic = dcic(1)/dtime
  ! dotcic = dotcic*conjg(cic(1))
  ! dotorb = 0d0
  ! call hprod_dotorb_tdcis(dtime, hwfn, h0orb, h1orb, dotorb)
  ! call hprod_mkovlp(xrad(nrad), orb, dotorb, o_dto)
  ! call hprod_mkovlp(xrad(nrad), orb, h0orb, o_h0)
  ! call hprod_mkovlp(xrad(nrad), orb, h1orb, o_h1)
  ! ! write(6, "('hprod_htot_tdcis:debug:orbitals')")
  ! ! do irad = 1, 40
  ! !    write(6, "(9es20.10)") xrad(irad), hwfn(irad, 0, 1), h0orb(irad, 1, 1), h1orb(irad, 0, 1), dotorb(irad, 0, 1)
  ! ! end do
  ! ! write(6 , "('hprod_htot_tdcis:debug', 10e20.10)") cic, o_o(1,1) dotcic, o_w(1,1), o_dto(1,1)
  ! write(6 , "('hprod_htot_tdcis:debug', 7e20.10)") real(dotcic), real(o_w(1,1)), real(o_dto(1,1)), o_w(1,1)/dtime
  !debug

! debug
! write(6, "('hprod_htot: wfn')")
! call hprod_printorb(wfn)
! write(6, "('hprod_htot: hwfn')")
! call hprod_printorb(hwfn)
! stop '@ hprod_htot for debug.'
! debug

end subroutine hprod_htot_tdcis
!#######################################################################
subroutine hprod_htotx_tdcis(dtime, lfield, wfn, cic)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : ecs_flag, irad_ecs
  use mod_control, only : exact3j
  use mod_bas, only : nbas, ngrid
  use mod_const, only : zero, czero, runit, cfalse, iunit
  use mod_hprod, only : rho1, rho2, v2sph, v2ang, v2ange, v2ango
  use mod_ormas, only : neltot, nfun, nfcore, ndcore, nact, lcic
  use mod_control, only : ioorot, isplit, iprojfc, igauge, dft_type, PSP
  use mod_hprod, only : ene_fcore, ene_dcore, ene_core, ene_act, ene_tot
  use mod_hprod, only : dcic, den1, den2, rden, rrden, den2r, int1e, int2e
  use mod_hprod, only : orb, orbg, h0orb, h1orb, gorb, gorbg, v2orb, v2orbg
  use mod_hprod, only : gorbe,gorbo,v2orbe,v2orbo,torb3j,rho23j,v2sph3j,orbe,orbo

  !tdcis
  use mod_sph, only : nlat
  use mod_rad, only : nrad
  use mod_control, only : tdcis_rvg
  use mod_hprod, only : orb0, orb0rot, h1orb0, orb0rotg, v2ang0
  use mod_hprod, only : ezphi
  use mod_ormas, only : nfcore_tdcis
  !tdcis

  !debug
  use mod_rad, only : xrad, nradfc
  !debug

  !debug
  !  use mod_hprod, only : projhigh
  !debug

  implicit none
  real(c_double), intent(in) :: dtime
  real(c_double), intent(in) :: lfield(1:3, 1:3)
  complex(c_double_complex), intent(in) :: wfn(1:*)
  complex(c_double_complex), intent(in) :: cic(1:lcic)

  integer(c_int) :: size_orb, size_forb, size_cic
  real(c_double), external :: hprod_ene_fcore
  real(c_double), external :: hprod_ene_dcore
  real(c_double), external :: hprod_ene_act
  complex(c_double_complex) :: zfield, efield, afield
  complex(c_double_complex) :: h1phichi(1:nfun, 1:nfun)
  complex(c_double_complex) :: ezphichi(1:nfun, 1:nfun)
  real(c_double) :: radmax

  ! debug
  integer(c_int) :: irad
  ! debug
  
  ! tdcis
  integer(c_int) :: ifun
  integer(c_int) :: fsize, osize, v2size
  fsize = nfun*nfun
  osize = nbas*nfun
  v2size = (nrad - 1)*nlat*nfun*nfun
  ! tdcis  

  ! laser field
  zfield = lfield(3, 1) ! chosen gauge
  efield = lfield(3, 2) ! electric field
  afield = lfield(3, 3) ! vector potential

  ! scratch initialization
  call hprod_htot_init

  ! current orbitals
  call hprod_orbin(lfield, wfn, orb, orbg)
  call hprod_orbin_tdcis(lfield)

  !debug
  ! call util_print_vec(osize, orb0, 'hprod.orb0')
  ! call util_print_vec(osize, orb0rot, 'hprod.orb0rot')
  ! call util_print_vec(v2size, v2ang0, 'hprod.v2ang0')
  ! stop
  !debug
!debug
!  stop 'hprod_htotx_tdcis: 1'
!debug

  ! RDMs
!Sato  call ormas_mkden1(cic, den1)
!Sato  call ormas_mkden2(cic, den1, den2)
!Sato  call hprod_invden(den1, rden, rrden)

  ! atomic hamiltonian
  if (ioorot == 0 .or. isplit == 0) then
     call hprod_tprod_dyn(orb, h0orb)
  end if
  if (iprojfc == 2) then
     call hprod_tprod_fc(orb, h0orb)
  end if

!debug
!  stop 'hprod_htotx_tdcis: 2'
!debug
  ! pseudopotential
  if (PSP) call hprod_projpp(runit, lfield, orb, h0orb)
!debug
!  stop 'hprod_htotx_tdcis: 3'
!debug

  ! laser hamiltonian
  if (ioorot .ne. 1 .or. isplit .ne. 1) then
    if (igauge == 0) then
       call hprod_zprod_dyn(zfield, orb, h1orb)
       call hprod_zprod_dyn(zfield, orb0rot, h1orb0)
    else
       if(tdcis_rvg) then
          call hprod_zprod_dyn(efield, orb0rot, ezphi)
       end if
       call hprod_pzprod_dyn(zfield, orb, h1orb)
       call hprod_pzprod_dyn(zfield, orb0rot, h1orb0) 
    end if
  end if

  if (iprojfc == 2) then
     if (igauge == 0) then
        call hprod_zprod_fc(zfield, orb, h1orb)
        call hprod_zprod_fc(zfield, orb0rot, h1orb0)
     else
        if(tdcis_rvg) then
          call hprod_zprod_fc(efield, orb0rot, ezphi)
        end if
        call hprod_pzprod_fc(zfield, orb, h1orb)
        call hprod_pzprod_fc(zfield, orb0rot, h1orb0)
     end if
  end if
  !NEW
  if (nfcore_tdcis > 0 .and. igauge == 1) call hprod_zprod_fc(efield, orb, h1orb);
  ! if (nfcore > 0 .and. igauge == 1) call hprod_zprod_fc(efield, orb, h1orb);
  !NEW

  ! debug
  ! if(real(zfield) .ne. 0d0) then
  !    write(6,'("###hprod_htot_tdcis:debug### zfield = ", 2f20.10)') zfield
  !    do irad = 1, nradfc
  !       ! write(6,"(7f20.10)") xrad(irad), orb0rot(irad,0,1), h1orb0(irad,0,1), h1orb0(irad,1,1)
  !       write(6,'(f10.5, 4es20.10)') xrad(irad), orb0rot(irad,0,1), h1orb0(irad,1,1)
  !    end do
  ! end if
  ! debug

  ! meanfield operator
  if (neltot(3) >= 2) then
     ! Orimo_ECS
!tdcis_sato     if (ecs_flag == 0) then
     !call hprod_mkv2mf2_dyn
     call hprod_mkrho2_tdcis
     call hprod_mkv2mf_tdcis
!tdcis_sato     else
!tdcis_sato        write(6, "('ecs_option: nyi')")
!tdcis_sato        stop 
!tdcis_sato     end if
     ! Orimo_ECS
     !call hprod_mfprod
     call hprod_mfprod2_tdcis(orb0rotg, orbg, v2ang, v2ang0, gorbg)
     call bas_ang2sph1_dyn(gorbg, gorb)
  end if

  ! debug
  !  if (projhigh) then
  !     call hprod_projhigh(h0orb)
  !     call hprod_projhigh(h1orb)
  !     call hprod_projhigh(gorb)
  !  end if
  ! debug

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
!Sato  if (ndcore > 0) ene_dcore = hprod_ene_dcore(cfalse, orb, h0orb, h1orb, gorb)
!Sato  if (nact > 0) ene_act = hprod_ene_act(int1e, int2e, den1, den2)
!Sato  ene_core = ene_fcore + ene_dcore
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

!tdcis_sato
  if (ecs_flag == 0) then
     radmax = xrad(nrad)
  else 
     radmax = xrad(irad_ecs-1)
  end if
!tdcis_sato

  ! ci coefficient
!tdcis_sato  if (ecs_flag == 0) then
  dcic(1) = 0d0
  if(.not. tdcis_rvg) then 
     !  LG or VG
     call hprod_mkovlp(radmax, h1orb0, orb, h1phichi) !<phi|h1|chi>
     ! do ifun = nfcore + 1, nfun
     do ifun = nfcore_tdcis + 1, nfun
        dcic(1) = dcic(1) + h1phichi(ifun,ifun)
     end do
  else 
     ! rVG
     call hprod_mkovlp(radmax, ezphi, orb, ezphichi) !<phi|Er|phi>
     ! do ifun = nfcore + 1, nfun
     do ifun = nfcore_tdcis + 1, nfun
        dcic(1) = dcic(1) + ezphichi(ifun,ifun)
     end do
  end if
  dcic(1) = dcic(1)*(-iunit)*sqrt(2d0)*dtime
!tdcis_sato  else
!tdcis_sato     write(6, "('ecs_option: not yet implemented')")
!tdcis_sato     stop 
!tdcis_sato     call hprod_htot_mkxmat_fc_ecs(dtime, orb, h0orb, h1orb, gorb, v2orb)
!tdcis_sato     call hprod_htot_mkxmat_cc_ecs(dtime, orb, h0orb, h1orb, gorb, v2orb)
!tdcis_sato     call hprod_htot_mkxmat_ca_ecs(dtime, orb, h0orb, h1orb, gorb, v2orb)
!tdcis_sato     call hprod_htot_mkxmat_aa1_ecs(dtime, orb, h0orb, h1orb, gorb, v2orb)
!tdcis_sato     call hprod_htot_mkxmat_aa2_ecs(dtime, orb, h0orb, h1orb, gorb, v2orb, cic, dcic)
!tdcis_sato  end if

end subroutine hprod_htotx_tdcis
!#######################################################################
