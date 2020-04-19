!###################################################################################
integer(c_int) function tdcc_spin_fock(s1, s2)

  use, intrinsic :: iso_c_binding

!  1 : f(p1,q1)
!  2 : f(p2,q2)
!  0 : others

  implicit none
  integer(c_int), intent(in) :: s1, s2
  integer(c_int) :: spin_type

  if (s1 .ne. s2) then
     spin_type = 0
  else if (s1 == 1 .and. s2 == 1) then
     spin_type = 1
  else if (s1 == 2 .and. s2 == 2) then
     spin_type = 2
  else
     spin_type = 0
  end if

  tdcc_spin_fock = spin_type

end function tdcc_spin_fock
!###################################################################################
