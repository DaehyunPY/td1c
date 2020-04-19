!######################################################################
subroutine ormas_hcic1(int1e, cic, hcic)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nact, cic_old

  implicit none
  complex(c_double_complex), intent(in) :: int1e(1:nact, 1:nact)
  complex(c_double_complex), intent(in) :: cic(1:*)
  complex(c_double_complex), intent(inout) :: hcic(1:*)

  if (.not. cic_old) then
     call ormas_hcic1_ras(int1e, cic, hcic)
  else
     call ormas_hcic1_old(int1e, cic, hcic)
  end if

end subroutine ormas_hcic1
!######################################################################
subroutine ormas_hcic1p(eref, int1e, cic, hcic)

  use, intrinsic :: iso_c_binding
  use mod_const, only : czero, runit
  use mod_ormas, only : nact, lcic, ormas_phase, cic_old

  implicit none
  real(c_double), intent(in) :: eref
  complex(c_double_complex), intent(in) :: int1e(1:nact, 1:nact)
  complex(c_double_complex), intent(in) :: cic(1:*)
  complex(c_double_complex), intent(inout) :: hcic(1:*)
  complex(c_double_complex) :: phase

  phase = -runit * eref
  call ormas_hcic1(int1e, cic, hcic)
  if (ormas_phase) then
     call zaxpy_omp(lcic, phase, cic, hcic)
  end if

end subroutine ormas_hcic1p
!######################################################################
subroutine ormas_zhcic1(zfac, int1e, cic, hcic)

  use, intrinsic :: iso_c_binding
  use mod_const, only : czero
  use mod_ormas, only : nact, lcic, cic_old

  implicit none
  complex(c_double_complex), intent(in) :: zfac
  complex(c_double_complex), intent(in) :: int1e(1:nact, 1:nact)
  complex(c_double_complex), intent(in) :: cic(1:*)
  complex(c_double_complex), intent(inout) :: hcic(1:*)
  complex(c_double_complex), allocatable :: tcic(:)

  allocate(tcic(1:lcic))
  tcic = czero
  call ormas_hcic1(int1e, cic, tcic)
  call zaxpy_omp(lcic, zfac, tcic, hcic)

  deallocate(tcic)

end subroutine ormas_zhcic1
!######################################################################
