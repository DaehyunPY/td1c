!######################################################################
subroutine ormas_spin(cic, sz, s2)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : cic_old

  implicit none
  complex(c_double_complex), intent(in) :: cic(1:*)
  real(c_double), intent(out) :: sz, s2

  if (.not. cic_old) then
     call ormas_spin_ras(cic, sz, s2)
  else
     call ormas_spin_old(cic, sz, s2)
  end if

end subroutine ormas_spin
!######################################################################
