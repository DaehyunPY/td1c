!######################################################################
subroutine ormas_mkden2(cic, den1, den2)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : smul
  use mod_ormas, only : den2_type,nact,nelact,fab_den2,cic_old,tdcc

  implicit none
  complex(c_double_complex) , intent(in) :: cic(1:*)
  complex(c_double_complex) , intent(in) :: den1(1:*)
  complex(c_double_complex) , intent(out) :: den2(1:*)

  if (nact == 0) return

  if (tdcc) then
     call tdcc_mkden2(cic, den2)
  else if (den2_type == 0) then
     if (.not. cic_old) then
        call ormas_mkden2_ras(cic, den1, den2)
     else
        call ormas_mkden2_old(cic, den1, den2)
     end if
!  else if (den2_type == 1) then
!     call ormas_mkden2_v1(cic, den1, den2)
!  else if (den2_type == 2) then ! previous default
!     call ormas_mkden2_v2(cic, den1, den2)
!  else if (den2_type == 3) then
!     call ormas_mkden2_v3(cic, den1, den2)
!  else if (den2_type == 4) then
!     if (fab_den2 .and. smul == 1 .and. nelact(1) == nelact(2)) then
!        call ormas_mkden2_v4(cic, den1, den2) ! new version by fabian
!     else
!        call ormas_mkden2_v2(cic, den1, den2)
!     end if
  else
     write(6, "('bad den2_type.')")
     stop
  end if

end subroutine ormas_mkden2
!######################################################################
subroutine ormas_print_den2(den2)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : smul
  use mod_ormas, only : den2_type, nact, nelact, fab_den2, cic_old

  implicit none
  complex(c_double_complex) , intent(in) :: den2(1:nact,1:nact,1:nact,1:nact)
  integer(c_long) :: i,j,k,l

  if (nact == 0) return
  do i = 1, nact
  do j = 1, nact
  do k = 1, nact
  do l = 1, nact
     write(6,"(4i5,2f10.5)") i,j,k,l,den2(i,j,k,l)
  end do
  end do
  end do
  end do

end subroutine ormas_print_den2
!######################################################################
