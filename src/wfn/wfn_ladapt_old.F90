!///////////////////////////////////////////////////////////////////////
subroutine wfn_ladapt_old(orb)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad
  use mod_sph, only : lmax1
  use mod_bas, only : lval
  use mod_ormas, only : nfun
  use mod_const, only : czero

  implicit none
  complex(c_double_complex), intent(inout) :: orb(1:(nrad-1), 0:lmax1, 1:nfun)

  integer(c_int) :: ifun, l, li, irad

  do ifun = 1, nfun
     li = lval(ifun)
     do l = 0, lmax1
        if (l .ne. li) then
           do irad = 1, nrad - 1
              orb(irad, l, ifun) = czero
           end do
        end if
     end do
  end do

end subroutine wfn_ladapt_old
!///////////////////////////////////////////////////////////////////////
subroutine wfn_ladapt_old_core(orb)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad
  use mod_sph, only : lmax1
  use mod_bas, only : lval
  use mod_ormas, only : ncore, nfun
  use mod_const, only : czero

  implicit none
  complex(c_double_complex), intent(inout) :: orb(1:(nrad-1), 0:lmax1, 1:nfun)

  integer(c_int) :: ifun, l, li, irad

  do ifun = 1, ncore
     li = lval(ifun)
     do l = 0, lmax1
        if (l .ne. li) then
           do irad = 1, nrad - 1
              orb(irad, l, ifun) = czero
           end do
        end if
     end do
  end do

end subroutine wfn_ladapt_old_core
!///////////////////////////////////////////////////////////////////////
