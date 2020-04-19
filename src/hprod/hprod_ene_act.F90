!######################################################################
real(c_double) function hprod_ene_act(int1e, int2e, den1, den2)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nact
  use mod_const, only : half, czero

  implicit none
  complex(c_double_complex), intent(in) :: int1e(1:nact, 1:nact)
  complex(c_double_complex), intent(in) :: int2e(1:nact, 1:nact, 1:nact, 1:nact)
  complex(c_double_complex), intent(in) :: den1(1:nact, 1:nact)
  complex(c_double_complex), intent(in) :: den2(1:nact, 1:nact, 1:nact, 1:nact)

  integer(c_int) :: iact, jact, kact, lact
  complex(c_double_complex) :: tmp1, tmp2

!DEBUG
!  do iact = 1, nact
!     do jact = 1, nact
!        write(6, "('int1e: ', 2i5, 2f20.10)") iact, jact, int1e(iact, jact)
!     end do
!  end do
!  do iact = 1, nact
!     do jact = 1, nact
!        do kact = 1, nact
!           do lact = 1, nact
!              write(6, "('int2e: ', 4i5, 2f20.10)") iact, jact, kact, lact, int2e(iact, jact, kact, lact)
!           end do
!        end do
!     end do
!  end do
!DEBUG
  
  tmp1 = czero
  tmp2 = czero
  do iact = 1, nact
     do jact = 1, nact
        tmp1 = tmp1 + int1e(iact, jact) * den1 (jact, iact)
        do kact = 1, nact
           do lact = 1, nact
              tmp2 = tmp2 + int2e(iact, jact, kact, lact) &
                          * den2 (jact, iact, lact, kact)
           end do
        end do
     end do
  end do

!  write(6,"('hprod_ene_act: int1 = ', f20.10)") dble(int1e(1,1))
!  write(6,"('hprod_ene_act: den1 = ', f20.10)") dble(den1(1,1))
!  write(6,"('hprod_ene_act: ene1 = ', f20.10)") dble(tmp1)
  hprod_ene_act = dble(tmp1) + dble(tmp2) * half
  return

end function hprod_ene_act
!######################################################################
real(c_double) function hprod_ene_act1(int1e, den1)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nact
  use mod_const, only : half, czero

  implicit none
  complex(c_double_complex), intent(in) :: int1e(1:nact, 1:nact)
  complex(c_double_complex), intent(in) :: den1(1:nact, 1:nact)

  integer(c_int) :: iact, jact
  complex(c_double_complex) :: tmp1

  tmp1 = czero
  do iact = 1, nact
     do jact = 1, nact
        tmp1 = tmp1 + int1e(iact, jact) * den1 (jact, iact)
     end do
  end do

  hprod_ene_act1 = dble(tmp1)
  return

end function hprod_ene_act1
!######################################################################
