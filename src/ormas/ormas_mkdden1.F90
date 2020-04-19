!######################################################################
subroutine ormas_mkdden1(fac, cic, dcic, dden1)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : cic_old

  implicit none
  complex(c_double_complex), intent(in) :: fac
  complex(c_double_complex), intent(in) :: cic(1:*)
  complex(c_double_complex), intent(in) :: dcic(1:*)
  complex(c_double_complex), intent(inout) :: dden1(1:*)

  if (.not. cic_old) then
     call ormas_mkdden1_ras(fac, cic, dcic, dden1)
  else
     call ormas_mkdden1_old(fac, cic, dcic, dden1)     
  end if

end subroutine ormas_mkdden1
!######################################################################
