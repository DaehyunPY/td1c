!######################################################################
subroutine ormas_mkden1(cic, den1)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : den1_type,nact,cic_old,tdcc

  implicit none
  complex(c_double_complex), intent(in) :: cic(1:*)
  complex(c_double_complex), intent(out) :: den1(1:*)

  if (nact == 0) return

  if (tdcc) then
     call tdcc_mkden1(cic, den1)
  else if (den1_type == 0) then
     if (.not. cic_old) then
        call ormas_mkden1_ras(cic, den1)
     else
        call ormas_mkden1_old(cic, den1)
     end if
!  else if (den1_type == 1) then ! previous default
!     call ormas_mkden1_v1(cic, den1)
!  else if (den1_type == 2) then
!     call ormas_mkden1_v2(cic, den1)
  else
     write(6, "('bad den1_type.')")
     stop
  end if

end subroutine ormas_mkden1
!######################################################################
