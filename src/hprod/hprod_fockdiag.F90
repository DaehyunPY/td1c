!#######################################################################
subroutine hprod_fockdiag(lfield, wfn, cic, eig)

  use, intrinsic :: iso_c_binding
  use mod_const, only : one
  use mod_bas, only : nbas
  use mod_hprod, only : orb, h0orb, gorb
  use mod_ormas, only : nfun
!DEBUG
  use mod_rad, only : nrad
  use mod_sph, only : lmax1
!DEBUG

  implicit none
  real(c_double), intent(in) :: lfield(1:3, 1:3)
  complex(c_double_complex), intent(inout) :: wfn(1:*)
  complex(c_double_complex), intent(inout) :: cic(1:*)
  real(c_double), intent(out) :: eig(1:*)

!DEBUG
!  integer(c_long) :: ifun, l, irad
!DEBUG

  write(6,"('hprod_fockdiag: use gfock instead of mpfock')")
! call hprod_mpfock(lfield, wfn, cic)
  call hprod_gfock(lfield, wfn, cic)
! call hprod_htot(icaldt, dtime_xx, lfield, wfn, cic, hwfn_xx, hcic_xx)

!DEBUG
!  write(6, "('hprod_fockdiagx: orb')")
!  do ifun = 1, nfun
!     write(6, "('# ifun = ', i5)") ifun
!     do irad = 1, nrad - 1
!        do l = 0, lmax1
!           write(6, "(E15.5)", advance = 'no') dble(orb(irad, l, ifun))
!        end do
!        write(6, *)
!     end do
!     write(6, *)
!     write(6, *)
!  end do
!  write(6, "('hprod_fockdiagx: h0orb')")
!  do ifun = 1, nfun
!     write(6, "('# ifun = ', i5)") ifun
!     do irad = 1, nrad - 1
!        do l = 0, lmax1
!           write(6, "(E15.5)", advance = 'no') dble(h0orb(irad, l, ifun))
!        end do
!        write(6, *)
!     end do
!     write(6, *)
!     write(6, *)
!  end do
!DEBUG

  call hprod_fockdiagx(orb, gorb, eig)
  call util_zcopy(nfun*nbas, orb, 1, wfn, 1)

end subroutine hprod_fockdiag
!#######################################################################
subroutine hprod_fockdiagx(wfn, fwfn, eig)

  use, intrinsic :: iso_c_binding
  use mod_sph, only : lmax1
  use mod_const, only : czero
  use mod_ormas, only : nfun
  use mod_bas, only : nbas, nval, lval, mval

  implicit none
  complex(c_double_complex), intent(inout) :: wfn(1:nbas, 1:*)
  complex(c_double_complex), intent(inout) :: fwfn(1:nbas, 1:*)
  real(c_double), intent(out) :: eig(1:*)

  integer(c_long) :: ifun, jfun
  complex(c_double_complex), allocatable :: fock(:,:)
  complex(c_double_complex), allocatable :: umat(:,:)

  integer(c_long) :: l, num0, iorb, jorb
  integer(c_long) :: map_i(1:nfun), map_j(1:nfun), map_n(1:nfun)


  do l = 0, lmax1
     num0 = 0
     do ifun = 1, nfun
        if (lval(ifun) == l .and. mval(ifun) == 0) then
           num0 = num0 + 1
           map_i(num0) = ifun
           map_n(num0) = nval(ifun)
        end if
     end do

     if (num0 == 0) cycle

     do jfun = 1, nfun
        if (lval(jfun) == l .and. mval(jfun) .ne. 0) then
           do iorb = 1, num0
              ifun = map_i(iorb)
              if (nval(jfun) == map_n(iorb)) then
                 map_j(jfun) = ifun
              end if
           end do
        end if
     end do

     allocate(fock(1:num0, 1:num0))
     allocate(umat(1:num0, 1:num0))

     fock(1:num0, 1:num0) = czero
     do iorb = 1, num0
        ifun = map_i(iorb)
        fock(iorb, iorb) = dot_product(wfn(1:nbas, ifun), fwfn(1:nbas, ifun))
        do jorb = 1, iorb - 1
           jfun = map_i(jorb)
           fock(jorb, iorb) = dot_product(wfn(1:nbas, jfun), fwfn(1:nbas, ifun))
           fock(iorb, jorb) = conjg(fock(jorb, iorb))
        end do
     end do

!debug
!     write(6, "('fock matrix: l = ', i5)") l
!     do iorb = 1, num0
!        ifun = map_i(iorb)
!        do jorb = 1, num0
!           jfun = map_i(jorb)
!           write(6, "(f20.10)", advance = 'no') dble(fock(jorb, iorb))
!        end do
!        write(6, *)
!     end do
!debug

     call futil_diag_comp(.false., num0, fock, umat)

!debug
!     write(6, "('transformation matrix: l = ', i5)") l
!     do iorb = 1, num0
!        ifun = map_i(iorb)
!        do jorb = 1, num0
!           jfun = map_i(jorb)
!           write(6, "(f20.10)", advance = 'no') dble(umat(jorb, iorb))
!        end do
!        write(6, *)
!     end do
!debug

     do iorb = 1, num0
        ifun = map_i(iorb)
        eig(ifun) = dble(fock(iorb, iorb))
        fwfn(1:nbas, ifun) = wfn(1:nbas, ifun)
        wfn (1:nbas, ifun) = czero
     end do
     do jfun = 1, nfun
        if (lval(jfun) == l .and. mval(jfun) .ne. 0) then
           eig(jfun) = eig(map_j(jfun))
           wfn(1:nbas, jfun) = czero
        end if
     end do

     do iorb = 1, num0
        ifun = map_i(iorb)
        do jorb = 1, num0
           jfun = map_i(jorb)
           wfn(1:nbas, ifun) = wfn(1:nbas, ifun) + fwfn(1:nbas, jfun) * umat(jorb, iorb)
        end do
     end do

     do jfun = 1, nfun
        if (lval(jfun) == l .and. mval(jfun) .ne. 0) then
           wfn(1:nbas, jfun) = wfn(1:nbas, map_j(jfun))
        end if
     end do

     deallocate(umat)
     deallocate(fock)
  end do

end subroutine hprod_fockdiagx
!#######################################################################
