!#######################################################################
subroutine hprod_gfock(lfield, wfn, cic)

  use, intrinsic :: iso_c_binding
  use mod_const, only : runit
  use mod_ormas, only : nfcore, neltot
  use mod_control, only : igauge, dft_type, PSP, exact3j
  use mod_hprod, only : orb, orbg, h0orb, h1orb, gorb, gorbg, v2orb, v2orbg
  use mod_hprod, only : rho1, rho2, v2sph, v2ang, den1, den2, rden, rrden
  use mod_rad, only : ecs_flag

  implicit none
  real(c_double), intent(in) :: lfield(1:3, 1:3)
  complex(c_double_complex), intent(in) :: wfn(*)
  complex(c_double_complex), intent(in) :: cic(*)
  complex(c_double_complex) :: zfield

  call hprod_htot_init                     ! scratch initialization
  call ormas_mkden1(cic, den1)             ! 1RDM
  call ormas_mkden2(cic, den1, den2)       ! 2RDM
  call hprod_invden(den1, rden, rrden)     ! 1RDM inverse
  call hprod_orbin(lfield, wfn, orb, orbg) ! current orbitals

  if (dft_type .ne. 0) then
     call hprod_gfock_dft(lfield, wfn, cic)
     return
  end if

  ! 1e operators
  zfield = lfield(3, 1)
! call hprod_tprod_ene(orb, h0orb)
  call hprod_tprod_all(orb, h0orb)
  if (PSP) call hprod_projpp(runit, lfield, orb, h0orb)

  if (igauge == 0) call hprod_zprod_all(zfield, orb, h1orb)
  if (igauge == 1) call hprod_pzprod_all(zfield, orb, h1orb)

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
        !call hprod_mfprod_gfock
        !call hprod_mfprod
        call hprod_mfprod2_gfock
        !call hprod_mfprod2
        call bas_ang2sph1_dyn(gorbg, gorb);
        call bas_ang2sph1_dyn(v2orbg, v2orb);
        if (nfcore > 0) then
           call bas_ang2sph1_fc(gorbg, gorb)
           call bas_ang2sph1_fc(v2orbg, v2orb)
        end if
     else
        call hprod_mkrho2_x3j_dyn
        call hprod_mkv2mf_x3j_dyn
        call hprod_mfprodx3j_gfock
     end if
  end if

  ! gfock vector
  call hprod_gfock_sum()
!  if (projhigh) then
!     call hprod_projhigh(gorb)
!  end if

!DEBUG
!  use mod_rad, only : nrad
!  use mod_sph, only : lmax1, lmax2
!  integer(c_int) :: isph, irad, ifun
!  write(6, "('gorb and v2orb:')")
!  do isph = 0, lmax1
!     do irad = 1, nrad - 1
!        write(6, "(2i5)", advance = 'no') isph, irad
!        do ifun = 1, nfun
!           write(6, "(4f20.10)", advance = 'no') gorb(irad, isph, ifun), v2orb(irad, isph, ifun)
!        end do
!        write(6, *)
!     end do
!  end do
!  write(6, "('rho2:')")
!  do isph = 0, lmax2
!     do irad = 1, nrad - 1
!        write(6, "(2i5)", advance = 'no') isph, irad
!        write(6, "(8f20.10)", advance = 'no') &
!  rho2(irad, isph, 1,1), rho2(irad, isph, 1,2), rho2(irad, isph, 2,1), rho2(irad, isph, 2,2)
!        write(6, *)
!     end do
!  end do
!  write(6, "('v2sph:')")
!  do isph = 0, lmax2
!     do irad = 1, nrad - 1
!        write(6, "(2i5)", advance = 'no') isph, irad
!        write(6, "(8f20.10)", advance = 'no') &
! v2sph(irad, isph, 1,1), v2sph(irad, isph, 1,2), v2sph(irad, isph, 2,1), v2sph(irad, isph, 2,2)
!        write(6, *)
!     end do
!  end do
!  stop
!DEBUG

end subroutine hprod_gfock
!#######################################################################
subroutine hprod_gfock_nofcx(lfield, wfn, cic)

  use, intrinsic :: iso_c_binding
  use mod_control, only : igauge
  use mod_ormas, only : nfcore, neltot
  use mod_hprod, only : orb, orbg, h0orb, h1orb, gorb, gorbg, v2orb, v2orbg
  use mod_hprod, only : rho1, rho2, v2sph, v2ang, den1, den2, rden, rrden
  use mod_rad, only : ecs_flag

  implicit none
  real(c_double), intent(in) :: lfield(1:3, 1:3)
  complex(c_double_complex), intent(in) :: wfn(*)
  complex(c_double_complex), intent(in) :: cic(*)
  complex(c_double_complex) :: zfield

  call hprod_htot_init                     ! scratch initialization
  call ormas_mkden1(cic, den1)             ! 1RDM
  call ormas_mkden2(cic, den1, den2)       ! 2RDM
  call hprod_invden(den1, rden, rrden)     ! 1RDM inverse
  call hprod_orbin(lfield, wfn, orb, orbg) ! current orbitals

  ! 1e operators
  zfield = lfield(3, 1)                    ! laser field
  call hprod_tprod_ene(orb, h0orb)         ! atomic hamiltonian
  if (igauge == 0) call hprod_zprod_all(zfield, orb, h1orb)  ! LG field
  if (igauge == 1) call hprod_pzprod_all(zfield, orb, h1orb) ! VG field

  ! meanfield operator
  if (neltot(3) >= 2) then
! Orimo_ECS
     if (ecs_flag == 0) then
        call hprod_mkrho2_dyn
        call hprod_mkv2mf_dyn
     else
        call hprod_mkrho2_dyn_ecs
        call hprod_mkv2mf_dyn_ecs
     end if
! Orimo_ECS
     call hprod_mfprod_gfock_nofcx
     call bas_ang2sph1_dyn(gorbg, gorb);
     call bas_ang2sph1_dyn(v2orbg, v2orb);
     if (nfcore > 0) then
        call bas_ang2sph1_fc(gorbg, gorb)
        call bas_ang2sph1_fc(v2orbg, v2orb)
     end if
  end if

  ! gfock vector
  call hprod_gfock_sum()

!DEBUG
!  use mod_rad, only : nrad
!  use mod_sph, only : lmax1, lmax2
!  integer(c_int) :: isph, irad, ifun
!  write(6, "('gorb and v2orb:')")
!  do isph = 0, lmax1
!     do irad = 1, nrad - 1
!        write(6, "(2i5)", advance = 'no') isph, irad
!        do ifun = 1, nfun
!           write(6, "(4f20.10)", advance = 'no') gorb(irad, isph, ifun), v2orb(irad, isph, ifun)
!        end do
!        write(6, *)
!     end do
!  end do
!  write(6, "('rho2:')")
!  do isph = 0, lmax2
!     do irad = 1, nrad - 1
!        write(6, "(2i5)", advance = 'no') isph, irad
!        write(6, "(8f20.10)", advance = 'no') &
!  rho2(irad, isph, 1,1), rho2(irad, isph, 1,2), rho2(irad, isph, 2,1), rho2(irad, isph, 2,2)
!        write(6, *)
!     end do
!  end do
!  write(6, "('v2sph:')")
!  do isph = 0, lmax2
!     do irad = 1, nrad - 1
!        write(6, "(2i5)", advance = 'no') isph, irad
!        write(6, "(8f20.10)", advance = 'no') &
!  v2sph(irad, isph, 1,1), v2sph(irad, isph, 1,2), v2sph(irad, isph, 2,1), v2sph(irad, isph, 2,2)
!        write(6, *)
!     end do
!  end do
!  stop
!DEBUG

end subroutine hprod_gfock_nofcx
!######################################################################
subroutine hprod_gfock_sum()

  use, intrinsic :: iso_c_binding
  use mod_sph, only : lmax1
  use mod_rad, only : nrad, nradfc
  use mod_ormas, only : nfcore, nfun
  use mod_hprod, only : h0orb, h1orb, gorb, v2orb

  implicit none
  integer(c_int) :: ifun, irad, l, llrf, ulrf, llrd, ulrd

  !$omp parallel default(shared) private(llrf, ulrf, llrd, ulrd)
  !###########################
  call util_omp_disp(1, nradfc,   llrf, ulrf)
  call util_omp_disp(1, nrad - 1, llrd, ulrd)
  do ifun = 1, nfcore
     do l = 0, lmax1
        do irad = llrf, ulrf
           gorb(irad, l, ifun) = v2orb(irad, l, ifun) + gorb(irad, l, ifun) + h1orb(irad, l, ifun) + h0orb(irad, l, ifun)
        end do
     end do
  end do
  do ifun = nfcore + 1, nfun
     do l = 0, lmax1
        do irad = llrd, ulrd
           gorb(irad, l, ifun) = v2orb(irad, l, ifun) + gorb(irad, l, ifun) + h1orb(irad, l, ifun) + h0orb(irad, l, ifun)
        end do
     end do
  end do
  !###########################
  !$omp end parallel

end subroutine hprod_gfock_sum
!#######################################################################
