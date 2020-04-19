!#######################################################################
subroutine hprod_orbene(lfield, wfn, cic, eig)

  use, intrinsic :: iso_c_binding
  use mod_const, only : one
  use mod_bas, only : nbas
  use mod_hprod, only : orb, gorb
  use mod_ormas, only : nfun

  implicit none
  real(c_double), intent(in) :: lfield(1:3, 1:3)
  complex(c_double_complex), intent(in) :: wfn(1:*)
  complex(c_double_complex), intent(in) :: cic(1:*)
  real(c_double), intent(out) :: eig(1:*)
  real(c_double), external :: hprod_ene_dft1

  write(6,"('hprod_orbene: use gfock instead of mpfock')")
! call hprod_mpfock(lfield, wfn, cic)
  call hprod_gfock(lfield, wfn, cic)

!debug
!  write(6, "('hprod_orbene:', f20.10)") hprod_ene_dft1(orb, gorb) * 0.5d+0
!debug
  call hprod_orbenex(.false., orb, gorb, eig)

end subroutine hprod_orbene
!#######################################################################
subroutine hprod_orbenex(fock_diag, wfn, fwfn, eig)

  use, intrinsic :: iso_c_binding
  use mod_const, only : czero
  use mod_bas, only : nbas, mval
  use mod_hprod, only : ovlp, xmat
  use mod_ormas, only : nfun

  implicit none
  logical, intent(in) :: fock_diag
  complex(c_double_complex), intent(inout) :: wfn(1:nbas, 1:*)
  complex(c_double_complex), intent(inout) :: fwfn(1:nbas, 1:*)
  real(c_double), intent(out) :: eig(1:*)

  integer(c_int) :: ifun, jfun
  real(c_double), external :: hprod_ene_dft1
  complex(c_double_complex), external :: util_zdotc

!debug
!  integer(c_int) :: ibas
!  write(6, "('hprod_orbenex: wfn and fwfn.')")
!  do ifun = 1, nfun
!     do ibas = 1, nbas
!        write(6, "(2i5, 4e20.10)") ifun, ibas, wfn(ibas, ifun), fwfn(ibas, ifun)
!     end do
!  end do
!  write(6, "('hprod_orbenex: wfn')")
!  call hprod_printorb(wfn)
!  write(6, "('hprod_orbenex: fwfn')")
!  call hprod_printorb(fwfn)
!debug

!debug
!  write(6, "('hprod_orbenex:', f20.10)") hprod_ene_dft1(wfn, fwfn) * 0.5d+0
!debug

  ovlp(1:nfun, 1:nfun) = czero
  do ifun = 1, nfun
!    ovlp(ifun, ifun) = dot_product(wfn(1:nbas, ifun), fwfn(1:nbas, ifun))
     ovlp(ifun, ifun) = util_zdotc(nbas, wfn(1, ifun), 1, fwfn(1, ifun), 1)
     do jfun = 1, ifun - 1
        if (mval(ifun) == mval(jfun)) then
!          ovlp(jfun, ifun) = dot_product(wfn(1:nbas, jfun), fwfn(1:nbas, ifun))
           ovlp(jfun, ifun) = util_zdotc(nbas, wfn(1, jfun), 1, fwfn(1, ifun), 1)
           ovlp(ifun, jfun) = conjg(ovlp(jfun, ifun))
        end if
     end do
  end do

!debug
  write(6, "('fock matrix:')")
  do ifun = 1, nfun
     do jfun = 1, nfun
        write(6, "(e20.10)", advance = 'no') dble(ovlp(jfun, ifun))
     end do
     write(6, *)
  end do
!debug

  call futil_diag_comp(.false., nfun, ovlp, xmat)

!debug
!  write(6, "('transformation matrix:')")
!  do ifun = 1, nfun
!     do jfun = 1, nfun
!        write(6, "(f20.10)", advance = 'no') dble(xmat(jfun, ifun))
!     end do
!     write(6, *)
!  end do
!debug

  do ifun = 1, nfun
     eig(ifun) = dble(ovlp(ifun, ifun))
  end do

  if (fock_diag) then
     fwfn(1:nbas, 1:nfun) = wfn(1:nbas, 1:nfun)
     wfn (1:nbas, 1:nfun) = czero
     do ifun = 1, nfun
        do jfun = 1, nfun
           wfn(1:nbas, ifun) = wfn(1:nbas, ifun) + fwfn(1:nbas, jfun) * xmat(jfun, ifun)
        end do
     end do
  end if

end subroutine hprod_orbenex
!#######################################################################
