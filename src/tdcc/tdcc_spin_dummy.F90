!###################################################################################
integer function tdcc_spin_dummy1(s1, s2)

  implicit none
  integer, intent(in) :: s1, s2
  integer, external :: tdcc_spin_fock

  tdcc_spin_dummy1 = 1
!  tdcc_spin_dummy1 = tdcc_spin_fock(s1, s2)

end function tdcc_spin_dummy1
!###################################################################################
integer function tdcc_spin_dummy2(s1, s2, s3, s4)

  implicit none
  integer, intent(in) :: s1, s2, s3, s4
  integer, external :: tdcc_spin_int2x

  tdcc_spin_dummy2 = 1
!  tdcc_spin_dummy2 = tdcc_spin_int2x(s1, s2, s3, s4)

end function tdcc_spin_dummy2
!###################################################################################
integer function tdcc_spin_dummy3(s1, s2, s3, s4, s5, s6)

  use mod_ormas, only : nelact

  implicit none
  integer, intent(in) :: s1, s2, s3, s4, s5, s6
  integer :: spin_type

!  stop 'tdcc_spin_dummy3: NEVER call me!'
  if (nelact(3) < 3) then
     spin_type = 0
  else if (s1==1.and.s2==1.and.s3==1.and.s4==1.and.s5==1.and.s6==1) then
     if (nelact(1) < 3) then
        spin_type = 0
     else
        spin_type = 1
     end if
  else if (s1==2.and.s2==2.and.s3==2.and.s4==2.and.s5==2.and.s6==2) then
     if (nelact(2) < 3) then
        spin_type = 0
     else
        spin_type = 1
     end if
  else
     spin_type = 1
  end if

  tdcc_spin_dummy3 = spin_type

end function tdcc_spin_dummy3
!###################################################################################
