!###################################################################################
integer(c_int) function tdcc_spin_t1inp(sp2, sh1)

!  1 : t(a1,i1)
!  2 : t(a2,i2)
!  0 : others

  use, intrinsic :: iso_c_binding

  implicit none
  integer(c_int), intent(in) :: sh1, sp2
  integer(c_int) :: spin_type

  if (sh1 .ne. sp2) then
     spin_type = 0
  else if (sh1 == 1 .and. sp2 == 1) then
     spin_type = 1
  else if (sh1 == 2 .and. sp2 == 2) then
     spin_type = 2
  else
     spin_type = 0
  end if

  tdcc_spin_t1inp = spin_type

end function tdcc_spin_t1inp
!###################################################################################
