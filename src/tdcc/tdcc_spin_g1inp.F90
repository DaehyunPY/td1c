!###################################################################################
integer function tdcc_spin_g1inp(sh1, sp2)

!  1 : g(i1,a1)
!  2 : g(i2,a2)
!  0 : others

  implicit none
  integer, intent(in) :: sh1, sp2
  integer :: spin_type

  if (sh1 .ne. sp2) then
     spin_type = 0
  else if (sh1 == 1 .and. sp2 == 1) then
     spin_type = 1
  else if (sh1 == 2 .and. sp2 == 2) then
     spin_type = 2
  else
     spin_type = 0
  end if

  tdcc_spin_g1inp = spin_type

end function tdcc_spin_g1inp
!###################################################################################
