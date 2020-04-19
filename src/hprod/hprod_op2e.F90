!######################################################################
subroutine hprod_op2e(rfac, int2e, den2, op2e)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : ncore, nact
  use mod_bas, only : mval

  implicit none
  real(c_double), intent(in) :: rfac
  complex(c_double_complex), intent(in) :: int2e(1:nact, 1:nact, 1:nact, 1:nact)
  complex(c_double_complex), intent(in) :: den2(1:nact, 1:nact, 1:nact, 1:nact)
  real(c_double), intent(out) :: op2e

  integer(c_int) :: iact, jact, kact, lact
  integer(c_int) :: mi, mj, mk, ml
  complex(c_double_complex) :: tmp
  complex(c_double_complex), parameter :: czero = (0.d+0, 0.d+0)
  
  tmp = czero
  do iact = 1, nact
     mi = mval(ncore + iact)
     do jact = 1, nact
        mj = mval(ncore + jact)
        do kact = 1, nact
           mk = mval(ncore + kact)
           do lact = 1, nact
              ml = mval(ncore + lact)
              if (mi + mk == mj + ml) then
                 tmp = tmp + int2e(iact, jact, kact, lact) &
                           * den2 (jact, iact, lact, kact)
              end if
           end do
        end do
     end do
  end do

  op2e = rfac * dble(tmp)

  !DEBUG
!  write(6, "('hprod_op2e. den2 (R)')")
!  do iact = 3, nact
!  do jact = 3, nact
!     do kact = 3, nact
!     do lact = 3, nact
!        write(6, "(f10.5)", advance = 'no') dble(den2(lact, kact, jact, iact))
!     end do
!     end do
!     write(6, *)
!  end do
!  end do
!!  write(6, "('hprod_op2e. den2 (I)')")
!!  do iact = 3, nact
!!  do jact = 3, nact
!!     do kact = 3, nact
!!     do lact = 3, nact
!!        write(6, "(f10.5)", advance = 'no') aimag(den2(lact, kact, jact, iact))
!!     end do
!!     end do
!!     write(6, *)
!!  end do
!!  end do
!
!  write(6, "('hprod_op2e. int2e (R)')")
!  do iact = 3, nact
!  do jact = 3, nact
!     do kact = 3, nact
!     do lact = 3, nact
!        write(6, "(f10.5)", advance = 'no') dble(int2e(kact, lact, iact, jact))
!     end do
!     end do
!     write(6, *)
!  end do
!  end do
!!  write(6, "('hprod_op2e. int2e (I)')")
!!  do iact = 3, nact
!!  do jact = 3, nact
!!     do kact = 3, nact
!!     do lact = 3, nact
!!        write(6, "(f10.5)", advance = 'no') aimag(int2e(kact, lact, iact, jact))
!!     end do
!!     end do
!!     write(6, *)
!!  end do
!!  end do
!
!  write(6, "('hprod_op2e. den2 (R) * int2e (R)')")
!  do iact = 3, nact
!  do jact = 3, nact
!     do kact = 3, nact
!     do lact = 3, nact
!        write(6, "(f10.5)", advance = 'no') dble(den2(lact, kact, jact, iact) * int2e(kact, lact, iact, jact))
!     end do
!     end do
!     write(6, *)
!  end do
!  end do
  !DEBUG

end subroutine hprod_op2e
!######################################################################
