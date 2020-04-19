!######################################################################
subroutine ormas_hcic(int1e, int2e, cic, hcic, ene_act)

  use, intrinsic :: iso_c_binding
  use mod_const, only : czero, runit
  use mod_ormas, only : cic_old, tdcc, nact, hcic_type, lcic, ormas_phase, ormas_donly
  use mod_cc, only : fock,int2x

  implicit none
  complex(c_double_complex), intent(in) :: int1e(1:*)
  complex(c_double_complex), intent(in) :: int2e(1:*)
  complex(c_double_complex), intent(in) :: cic(1:*)
  complex(c_double_complex), intent(inout) :: hcic(1:*)
  real(c_double), intent(in) :: ene_act
  complex(c_double_complex) :: phase

  if (nact == 0) return

  if (tdcc) then
     call tdcc_hcc12(int1e, int2e, cic, hcic)
     return
  else if (hcic_type == 0) then
     if (.not. cic_old) then
        call ormas_hcic_ras(int1e, int2e, cic, hcic)
!        call ormas_hcic_rasmat(int1e, int2e, cic, hcic)
     else
        call ormas_hcic_old(int1e, int2e, cic, hcic)
     end if
!     call ormas_hcic_v2(int1e, int2e, cic, hcic)
!  else if (hcic_type == 1) then ! for cas
!     call ormas_hcic_v1(int1e, int2e, cic, hcic)
!  else if (hcic_type == 2) then ! for ras
!     call ormas_hcic_v2(int1e, int2e, cic, hcic)
!  else if (hcic_type == 3) then ! for general ormas, not fully tested
!     call ormas_hcic_v3(int1e, int2e, cic, hcic)
  else
     write(6, "('bad hcic_type.')")
     stop
  end if

  if (ormas_phase) then
     phase = -runit * ene_act
     call zaxpy_omp(lcic, phase, cic, hcic)
  end if

!DEBUG
!  call util_print_vec(nact**2, int1e, "test.int1e")
!  call util_print_vec(nact**4, int2e, "test.int2e")
!  call ormas_cic_print(cic, "test.ciin")
!  call ormas_cic_print(hcic, "test.hcic")
!  call ormas_cic_print_old(cic, "test.ciinold")
!  call ormas_cic_print_old(hcic, "test.hcicold")
!DEBUG

end subroutine ormas_hcic
!######################################################################
!######################################################################
subroutine ormas_zhcic(zfac, int1e, int2e, cic, hcic)

  use, intrinsic :: iso_c_binding
  use mod_const, only : zero, czero, runit
  use mod_ormas, only : nact, hcic_type, lcic

  implicit none
  complex(c_double_complex), intent(in) :: zfac
  complex(c_double_complex), intent(in) :: int1e(1:*)
  complex(c_double_complex), intent(in) :: int2e(1:*)
  complex(c_double_complex), intent(in) :: cic(1:*)
  complex(c_double_complex), intent(inout) :: hcic(1:*)
  complex(c_double_complex), allocatable :: tcic(:)

  allocate(tcic(1:lcic))
  tcic = czero

  call ormas_hcic(int1e, int2e, cic, tcic, zero)
  call zaxpy_omp(lcic, zfac, tcic, hcic)

  deallocate(tcic)

end subroutine ormas_zhcic
!######################################################################
subroutine ormas_zhcicp(zfac, eref, int1e, int2e, cic, hcic)

  use, intrinsic :: iso_c_binding
  use mod_const, only : czero, runit
  use mod_ormas, only : nact, hcic_type, lcic

  implicit none
  complex(c_double_complex), intent(in) :: zfac
  real(c_double), intent(in) :: eref
  complex(c_double_complex), intent(in) :: int1e(1:*)
  complex(c_double_complex), intent(in) :: int2e(1:*)
  complex(c_double_complex), intent(in) :: cic(1:*)
  complex(c_double_complex), intent(inout) :: hcic(1:*)
  complex(c_double_complex), allocatable :: tcic(:)

  allocate(tcic(1:lcic))
  tcic = czero

  call ormas_hcic(int1e, int2e, cic, tcic, eref)
  call zaxpy_omp(lcic, zfac, tcic, hcic)

  deallocate(tcic)

end subroutine ormas_zhcicp
!######################################################################
