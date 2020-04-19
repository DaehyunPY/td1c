!######################################################################
subroutine hprod_get_norb(den1, wfn, nwfn)

  use, intrinsic :: iso_c_binding
  use mod_const, only : czero
  use mod_rad, only : nrad
  use mod_sph, only : lmax1
  use mod_bas, only : mval
  use mod_ormas, only : ncore, nact, nocc, nfun

  implicit none
  complex(c_double_complex), intent(inout) :: den1(1:nact, 1:nact)
  complex(c_double_complex), intent(in) :: wfn(1:(nrad-1), 0:lmax1, 1:*)
  complex(c_double_complex), intent(out) :: nwfn(1:(nrad-1), 0:lmax1, 1:*)

  integer(c_int) :: iact, jact, ifun, jfun, irad, l, m, mi, mj
  complex(c_double_complex), allocatable :: umat(:,:)

  allocate(umat(1:nact, 1:nact))
  if (nact > 0) then
     umat(1:nact, 1:nact) = czero
     call lapack_zheev_dsc(nact, den1, umat)
  end if

  ! DEBUG
  write(6, "('hprod_get_norb: natural occupation number')")
  do ifun = 1, ncore
     write(6, "(i5, f20.10)") ifun, 2.0
  end do
  do iact = 1, nact
     ifun = ncore + iact
     write(6, "(i5, f20.10)") ifun, dble(den1(iact, iact))
  end do
  do ifun = nocc + 1, nfun
     write(6, "(i5, f20.10)") ifun, 0.0
  end do  

  !write(6, "('eigen-solutions of den1 (R):')")
  !do iact = 1, nact
  !   write(6, "(5x, i12)", advance = 'no') iact
  !end do
  !write(6, *)
  !do iact = 1, nact
  !   write(6, "(5x, f12.5)", advance = 'no') dble(den1(iact, iact))
  !end do
  !write(6, *)
  !do jact = 1, nact
  !   do iact = 1, nact
  !      write(6, "(i5, f12.5)", advance = 'no') jact, dble(umat(jact, iact))
  !   end do
  !   write(6, *)
  !end do
  !write(6, "('eigen-solutions of den1 (I):')")
  !do iact = 1, nact
  !   write(6, "(5x, i12)", advance = 'no') iact
  !end do
  !write(6, *)
  !do iact = 1, nact
  !   write(6, "(5x, f12.5)", advance = 'no') aimag(den1(iact, iact))
  !end do
  !write(6, *)
  !do jact = 1, nact
  !   do iact = 1, nact
  !      write(6, "(i5, f12.5)", advance = 'no') jact, aimag(umat(jact, iact))
  !   end do
  !   write(6, *)
  !end do
  ! DEBUG

  nwfn(1:(nrad-1), 0:lmax1, 1:nfun) = czero
  do ifun = 1, ncore
     do l = abs(mval(ifun)), lmax1
        do irad = 1, nrad - 1
           nwfn(irad, l, ifun) = wfn(irad, l, ifun)
        end do
     end do
  end do
  do iact = 1, nact
     ifun = ncore + iact
     mi = mval(ifun)

     do jact = 1, nact
        jfun = ncore + jact
        mj = mval(jfun)

!bug        if (mi == mj) then
!bug           m = mi
!bug
        m = 0
        do l = abs(m), lmax1
           do irad = 1, nrad - 1
              nwfn(irad, l, ifun) = nwfn(irad, l, ifun) &
                                   + wfn(irad, l, jfun) * umat(jact, iact)
           end do
        end do
!bug        end if
     end do
  end do
  do ifun = nocc + 1, nfun
     do l = abs(mval(ifun)), lmax1
        do irad = 1, nrad - 1
           nwfn(irad, l, ifun) = wfn(irad, l, ifun)
        end do
     end do
  end do

  deallocate(umat)

end subroutine hprod_get_norb
!######################################################################
!subroutine hprod_get_norb_msep(den1, wfn, nwfn, umat)
!
!  use, intrinsic :: iso_c_binding
!  use mod_const, only : czero
!  use mod_rad, only : nrad
!  use mod_sph, only : lmax1, mmax1
!  use mod_bas, only : mval
!  use mod_ormas, only : ncore, nact, nocc, nvir, nfun
!
!  implicit none
!  complex(c_double_complex), intent(inout) :: den1(1:nact, 1:nact)
!  complex(c_double_complex), intent(in) :: wfn(1:(nrad-1), 0:lmax1, 1:*)
!  complex(c_double_complex), intent(out) :: nwfn(1:(nrad-1), 0:lmax1, 1:*)
!  complex(c_double_complex), intent(out) :: umat(1:nact, 1:nact)
!
!  integer(c_int) :: iact, jact, ifun, jfun, irad, l, m, mi, mj
!  complex(c_double_complex) :: fac1, fac2
!  complex(c_int), allocatable :: map_msep(:)
!  complex(c_double_complex), allocatable :: umat_msep(:,:)
!
!  allocate(map_msep(1:nact))
!
!  do m = 0, mmax1
!     norb_msep = 0
!     do iact = 1, nact
!        ifun = ncore + iact
!        if (mval(ifun) == m) then
!           norb_msep = norb_msep + 1
!           map_msep(norb_msep) = ifun
!        end if
!     end do
!
!     allocate(den1_msep(1:norb_msep, 1:norb_msep))
!     allocate(umat_msep(1:norb_msep, 1:norb_msep))
!
!     do iorb = 1, norb_msep
!        iact = map_msep(iorb)
!        do jorb = 1, norb_msep
!           jact = map_msep(jorb)
!           den1_msep(jorb, iorb) = den1(jact, iact)
!        end do
!     end do
!
!     umat_msep(1:orb_msep, 1:norb_msep) = czero
!     call lapack_zheev_dsc(norb_msep, den1_msep, umat_msep)
!
!  end do
!
!  allocate(umat(1:nact, 1:nact))
!  if (nact > 0) then
!     umat(1:nact, 1:nact) = czero
!     call lapack_zheev_dsc(nact, den1, umat)
!  end if
!
!  ! DEBUG
!  write(6, "('hprod_get_norb: natural occupation number')")
!  do ifun = 1, ncore
!     write(6, "(i5, f20.10)") ifun, 2.0
!  end do
!  do iact = 1, nact
!     ifun = ncore + iact
!     write(6, "(i5, f20.10)") ifun, dble(den1(iact, iact))
!  end do
!  do ifun = nocc + 1, nfun
!     write(6, "(i5, f20.10)") ifun, 0.0
!  end do  
!
!  !write(6, "('eigen-solutions of den1 (R):')")
!  !do iact = 1, nact
!  !   write(6, "(5x, i12)", advance = 'no') iact
!  !end do
!  !write(6, *)
!  !do iact = 1, nact
!  !   write(6, "(5x, f12.5)", advance = 'no') dble(den1(iact, iact))
!  !end do
!  !write(6, *)
!  !do jact = 1, nact
!  !   do iact = 1, nact
!  !      write(6, "(i5, f12.5)", advance = 'no') jact, dble(umat(jact, iact))
!  !   end do
!  !   write(6, *)
!  !end do
!  !write(6, "('eigen-solutions of den1 (I):')")
!  !do iact = 1, nact
!  !   write(6, "(5x, i12)", advance = 'no') iact
!  !end do
!  !write(6, *)
!  !do iact = 1, nact
!  !   write(6, "(5x, f12.5)", advance = 'no') aimag(den1(iact, iact))
!  !end do
!  !write(6, *)
!  !do jact = 1, nact
!  !   do iact = 1, nact
!  !      write(6, "(i5, f12.5)", advance = 'no') jact, aimag(umat(jact, iact))
!  !   end do
!  !   write(6, *)
!  !end do
!  ! DEBUG
!
!  nwfn(1:(nrad-1), 0:lmax1, 1:nfun) = czero
!  do ifun = 1, ncore
!     do l = abs(mval(ifun)), lmax1
!        do irad = 1, nrad - 1
!           nwfn(irad, l, ifun) = wfn(irad, l, ifun)
!        end do
!     end do
!  end do
!  do iact = 1, nact
!     ifun = ncore + iact
!     mi = mval(ifun)
!
!     do jact = 1, nact
!        jfun = ncore + jact
!        mj = mval(jfun)
!
!!bug        if (mi == mj) then
!!bug           m = mi
!!bug
!        m = 0
!        do l = abs(m), lmax1
!           do irad = 1, nrad - 1
!              nwfn(irad, l, ifun) = nwfn(irad, l, ifun) &
!                                   + wfn(irad, l, jfun) * umat(jact, iact)
!           end do
!        end do
!!bug        end if
!     end do
!  end do
!  do ifun = nocc + 1, nfun
!     do l = abs(mval(ifun)), lmax1
!        do irad = 1, nrad - 1
!           nwfn(irad, l, ifun) = wfn(irad, l, ifun)
!        end do
!     end do
!  end do
!
!  deallocate(umat)
!
!end subroutine hprod_get_norb_msep
!######################################################################
