!###################################################################################
integer(c_int) function tdcc_spin_t3inp(sp4, sp5, sp6, sh1, sh2, sh3)

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
  integer(c_int), intent(in) :: sh1,sh2,sh3, sp4,sp5,sp6
  integer(c_int) :: spin_type

  if (nelact(3)<3 .or. sh1+sh2+sh3.ne.sp4+sp5+sp6) then
     spin_type = 0
  else if (sp4==1.and.sp5==1.and.sp6==1 .and. sh1==1.and.sh2==1.and.sh3==1) then
     ! 1 (+++|+++)
     if (nelact(1) < 3) then
        spin_type = 0
     else
        spin_type = 1
     end if
  else if (sp4==2.and.sp5==2.and.sp6==2 .and. sh1==2.and.sh2==2.and.sh3==2) then
     ! 2 (---|---)
     if (nelact(2) < 3) then
        spin_type = 0
     else
        spin_type = 2
     end if
  else if (sp4==1.and.sp5==1.and.sp6==2 .and. sh1==1.and.sh2==1.and.sh3==2) then
     ! 3 (++-|++-)
     spin_type = 3
  else if (sp4==1.and.sp5==1.and.sp6==2 .and. sh1==1.and.sh2==2.and.sh3==1) then
     ! 4 (++-|+-+)
     spin_type = 4
  else if (sp4==1.and.sp5==1.and.sp6==2 .and. sh1==2.and.sh2==1.and.sh3==1) then
     ! 5 (++-|-++)
     spin_type = 5

  else if (sp4==1.and.sp5==2.and.sp6==1 .and. sh1==1.and.sh2==1.and.sh3==2) then
     ! 6 (+-+|++-)
     spin_type = 6
  else if (sp4==1.and.sp5==2.and.sp6==1 .and. sh1==1.and.sh2==2.and.sh3==1) then
     ! 7 (+-+|+-+)
     spin_type = 7
  else if (sp4==1.and.sp5==2.and.sp6==1 .and. sh1==2.and.sh2==1.and.sh3==1) then
     ! 8 (+-+|-++)
     spin_type = 8

  else if (sp4==2.and.sp5==1.and.sp6==1 .and. sh1==1.and.sh2==1.and.sh3==2) then
     ! 9 (-++|++-)
     spin_type = 9
  else if (sp4==2.and.sp5==1.and.sp6==1 .and. sh1==1.and.sh2==2.and.sh3==1) then
     ! 10 (-++|+-+)
     spin_type = 10
  else if (sp4==2.and.sp5==1.and.sp6==1 .and. sh1==2.and.sh2==1.and.sh3==1) then
     ! 11 (-++|-++)
     spin_type = 11

  else if (sp4==2.and.sp5==2.and.sp6==1 .and. sh1==2.and.sh2==2.and.sh3==1) then
     ! 12 (--+|--+)
     spin_type = 12
  else if (sp4==2.and.sp5==2.and.sp6==1 .and. sh1==2.and.sh2==1.and.sh3==2) then
     ! 13 (--+|-+-)
     spin_type = 13
  else if (sp4==2.and.sp5==2.and.sp6==1 .and. sh1==1.and.sh2==2.and.sh3==2) then
     ! 14 (--+|+--)
     spin_type = 14

  else if (sp4==2.and.sp5==1.and.sp6==2 .and. sh1==2.and.sh2==2.and.sh3==1) then
     ! 15 (-+-|--+)
     spin_type = 15
  else if (sp4==2.and.sp5==1.and.sp6==2 .and. sh1==2.and.sh2==1.and.sh3==2) then
     ! 16 (-+-|-+-)
     spin_type = 16
  else if (sp4==2.and.sp5==1.and.sp6==2 .and. sh1==1.and.sh2==2.and.sh3==2) then
     ! 17 (-+-|+--)
     spin_type = 17

  else if (sp4==1.and.sp5==2.and.sp6==2 .and. sh1==2.and.sh2==2.and.sh3==1) then
     ! 18 (+--|--+)
     spin_type = 18
  else if (sp4==1.and.sp5==2.and.sp6==2 .and. sh1==2.and.sh2==1.and.sh3==2) then
     ! 19 (+--|-+-)
     spin_type = 19
  else if (sp4==1.and.sp5==2.and.sp6==2 .and. sh1==1.and.sh2==2.and.sh3==2) then
     ! 20 (+--|+--)
     spin_type = 20
  else
     spin_type = 0
  end if

  tdcc_spin_t3inp = spin_type

end function tdcc_spin_t3inp
!###################################################################################
