!######################################################################
subroutine ormas_mkden2x(cic, den2x)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : cic_old

  implicit none
  complex(c_double_complex), intent(in) :: cic(1:*)
  complex(c_double_complex), intent(out) :: den2x(1:*)

  if (.not. cic_old) then
     call ormas_mkden2x_ras(cic, den2x)
  else
     call ormas_mkden2x_old(cic, den2x)
  end if

end subroutine ormas_mkden2x
!######################################################################
