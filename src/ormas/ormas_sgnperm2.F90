!################################################################################
subroutine ormas_sgnperm2(i0,i1,i0x,i1x,sgn)
  use,intrinsic :: iso_c_binding
  use mod_ormas,only : ncore

  implicit none
  !--------------------------------------------------------------------
  integer(c_int),intent(in) :: i0,i1
  integer(c_int),intent(out) :: i0x,i1x,sgn
  !--------------------------------------------------------------------

  if (i0 < i1) then
     sgn = (-1)**(i1-1)
     i0x = i0
     i1x = i1
  else
     sgn = (-1)**i1
     i0x = i1
     i1x = i0
  end if

end subroutine ormas_sgnperm2
!################################################################################
subroutine ormas_sgnperm3(i0,i1,i2,i0x,i1x,i2x,sgn)
  use,intrinsic :: iso_c_binding
  use mod_ormas,only : ncore

  implicit none
  !--------------------------------------------------------------------
  integer(c_int),intent(in) :: i0,i1,i2
  integer(c_int),intent(out) :: i0x,i1x,i2x,sgn
  !--------------------------------------------------------------------

  if (i0 < i1) then
     i0x = i0
     i1x = i1
     i2x = i2
     sgn = (-1)**(i1+i2)
  else if (i0 < i2) then
     i0x = i1
     i1x = i0
     i2x = i2
     sgn = (-1)**(i1+i2-1)
  else
     i0x = i1
     i1x = i2
     i2x = i0
     sgn = (-1)**(i1+i2)
  end if

end subroutine ormas_sgnperm3
!################################################################################
subroutine ormas_sgnperm4(i0,i1,i2,i3,i0x,i1x,i2x,i3x,sgn)
  use,intrinsic :: iso_c_binding
  use mod_ormas,only : ncore

  implicit none
  !--------------------------------------------------------------------
  integer(c_int),intent(in) :: i0,i1,i2,i3
  integer(c_int),intent(out) :: i0x,i1x,i2x,i3x,sgn
  !--------------------------------------------------------------------

  if (i0 < i1) then
     i0x = i0
     i1x = i1
     i2x = i2
     i3x = i3
     sgn = (-1)**(i1+i2+i3-1)
  else if (i0 < i2) then
     i0x = i1
     i1x = i0
     i2x = i2
     i3x = i3
     sgn = (-1)**(i1+i2+i3)
  else if (i0 < i3) then
     i0x = i1
     i1x = i2
     i2x = i0
     i3x = i3
     sgn = (-1)**(i1+i2+i3-1)
  else
     i0x = i1
     i1x = i2
     i2x = i3
     i3x = i0
     sgn = (-1)**(i1+i2+i3)
  end if

end subroutine ormas_sgnperm4
!################################################################################
subroutine ormas_sgnperm5(i0,i1,i2,i3,i4,i0x,i1x,i2x,i3x,i4x,sgn)
  use,intrinsic :: iso_c_binding
  use mod_ormas,only : ncore

  implicit none
  !--------------------------------------------------------------------
  integer(c_int),intent(in) :: i0,i1,i2,i3,i4
  integer(c_int),intent(out) :: i0x,i1x,i2x,i3x,i4x,sgn
  !--------------------------------------------------------------------

  if (i0 < i1) then
     i0x = i0
     i1x = i1
     i2x = i2
     i3x = i3
     i4x = i4
     sgn = (-1)**(i1+i2+i3+i4)
  else if (i0 < i2) then
     i0x = i1
     i1x = i0
     i2x = i2
     i3x = i3
     i4x = i4
     sgn = (-1)**(i1+i2+i3+i4-1)
  else if (i0 < i3) then
     i0x = i1
     i1x = i2
     i2x = i0
     i3x = i3
     i4x = i4
     sgn = (-1)**(i1+i2+i3+i4)
  else if (i0 < i4) then
     i0x = i1
     i1x = i2
     i2x = i3
     i3x = i0
     i4x = i4
     sgn = (-1)**(i1+i2+i3+i4-1)
  else
     i0x = i1
     i1x = i2
     i2x = i3
     i3x = i4
     i4x = i0
     sgn = (-1)**(i1+i2+i3+i4)
  end if

end subroutine ormas_sgnperm5
!################################################################################
