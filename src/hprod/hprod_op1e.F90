!######################################################################
subroutine hprod_op1e(rfac, int1e, den1, op1e)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : ncore, nact
  use mod_bas, only : mval

  implicit none
  real(c_double), intent(in) :: rfac
  complex(c_double_complex), intent(in) :: int1e(1:nact, 1:nact)
  complex(c_double_complex), intent(in) :: den1(1:nact, 1:nact)
  real(c_double), intent(out) :: op1e

  integer(c_int) :: iact, jact
  integer(c_int) :: ifun, jfun
  complex(c_double_complex) :: tmp
  complex(c_double_complex), parameter :: czero = (0.d+0, 0.d+0)
  
  tmp = czero
  do iact = 1, nact
     ifun = ncore + iact
     do jact = 1, nact
        jfun = ncore + jact
        if (mval(jfun) == mval(ifun)) then
           tmp = tmp + int1e(iact, jact) * den1 (jact, iact)
        end if
     end do
  end do

  op1e = rfac * dble(tmp)

  !DEBUG
!  write(6, "('hprod_op1e. den1 (R)')")
!  do iact = 1, nact
!     do jact = 1, nact
!        write(6, "(f20.10)", advance = 'no') dble(den1(jact, iact))
!     end do
!     write(6, *)
!  end do
!  write(6, "('hprod_op1e. den1 (I)')")
!  do iact = 1, nact
!     do jact = 1, nact
!        write(6, "(f20.10)", advance = 'no') aimag(den1(jact, iact))
!     end do
!     write(6, *)
!  end do
!  write(6, "('hprod_op1e. int1e (R)')")
!  do iact = 1, nact
!     do jact = 1, nact
!        write(6, "(f20.10)", advance = 'no') dble(int1e(jact, iact))
!     end do
!     write(6, *)
!  end do
!  write(6, "('hprod_op1e. int1e (I)')")
!  do iact = 1, nact
!     do jact = 1, nact
!        write(6, "(f20.10)", advance = 'no') aimag(int1e(jact, iact))
!     end do
!     write(6, *)
!  end do
  !DEBUG

end subroutine hprod_op1e
!######################################################################
