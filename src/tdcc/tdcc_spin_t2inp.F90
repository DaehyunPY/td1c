!###################################################################################
integer(c_int) function tdcc_spin_t2inp(sp3, sp4, sh1, sh2)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nelact
!  1 : t(a1b1,i1j1)
!  2 : t(a2b2,i2j2)
!  3 : t(a1b2,i1j2)
!  4 : t(a2b1,i2j1)
!  5 : t(a1b2,i2j1)
!  6 : t(a2b1,i1j2)
!  0 : others

  implicit none
  integer(c_int), intent(in) :: sh1, sh2, sp3, sp4
  integer(c_int) :: spin_type

  if (nelact(3) < 2 .or. sh1 + sh2 .ne. sp3 + sp4) then
     spin_type = 0
  else if (sh1 == 1 .and. sh2 == 1 .and. sp3 == 1 .and. sp4 == 1) then
     ! 1 (p+p+|h+h+)
     if (nelact(1) < 2) then
        spin_type = 0
     else
        spin_type = 1
     end if
  else if (sh1 == 2 .and. sh2 == 2 .and. sp3 == 2 .and. sp4 == 2) then
     ! 2 (p-p-|h-h-)
     if (nelact(2) < 2) then
        spin_type = 0
     else
        spin_type = 2
     end if
  else if (sh1 == 1 .and. sh2 == 2 .and. sp3 == 1 .and. sp4 == 2) then
     ! 3 (p+p-|h+h-)
     spin_type = 3
  else if (sh1 == 2 .and. sh2 == 1 .and. sp3 == 2 .and. sp4 == 1) then
     ! 4 (p-p+|h-h+)
     spin_type = 4
  else if (sh1 == 2 .and. sh2 == 1 .and. sp3 == 1 .and. sp4 == 2) then
     ! 5 (p+p-|h-h+)
     spin_type = 5
  else if (sh1 == 1 .and. sh2 == 2 .and. sp3 == 2 .and. sp4 == 1) then
     ! 6 (p-p+|h+h-)
     spin_type = 6
  else
     spin_type = 0
  end if

  tdcc_spin_t2inp = spin_type

end function tdcc_spin_t2inp
!###################################################################################
