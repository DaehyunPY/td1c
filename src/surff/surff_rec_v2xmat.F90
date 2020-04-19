!######################################################################
subroutine surff_rec_v2xmat(v2xmat)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : mval
  use mod_ormas, only : nfun, nfcore
  use mod_hprod, only : xmat

  implicit none
  complex(c_double_complex), intent(out) :: v2xmat(1:nfun, 1:nfun)

  integer(c_int) :: ifun, jfun

  do ifun = nfcore + 1, nfun
     do jfun = 1, nfun
        if (mval(ifun) == mval(jfun)) then
           v2xmat(jfun, ifun) = xmat(jfun, ifun)
        end if
     end do
  end do

end subroutine surff_rec_v2xmat
!######################################################################

