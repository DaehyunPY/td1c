!######################################################################
subroutine hprod_norb(wfn, cic, nwfn)

  use, intrinsic :: iso_c_binding
  use mod_const, only : czero
  use mod_rad, only : nrad
  use mod_sph, only : lmax1
  use mod_bas, only : mval
  use mod_hprod, only : den1
  use mod_ormas, only : ncore, nact, nocc, nfun

  implicit none
  complex(c_double_complex), intent(in) :: wfn(1:(nrad-1), 0:lmax1, 1:*)
  complex(c_double_complex), intent(in) :: cic(1:*)
  complex(c_double_complex), intent(out) :: nwfn(1:(nrad-1), 0:lmax1, 1:*)

  integer(c_long) :: iact, jact, ifun, jfun, irad, l, m, mi, mj
  complex(c_double_complex), allocatable :: umat(:,:)
  real(c_double), allocatable :: occn(:)


!old  allocate(umat(1:nact, 1:nact))
  call ormas_mkden1(cic, den1)
  allocate(occn(1:nfun))
!DEBUG!old     call hprod_norb2(den1, umat, occn)
!DEBUG!old     call hprod_norb3(den1, umat, occn)
  call hprod_norb4(den1, wfn, nwfn, occn)
  deallocate(occn)

!old  if (nact > 0) then
!old     umat(1:nact, 1:nact) = czero
!old     call lapack_zheev_dsc(nact, den1, umat)
!old  end if
!old
!old  ! DEBUG
!old  write(6, "('hprod_norb: natural occupation number')")
!old  do ifun = 1, ncore
!old     write(6, "(i5, e25.15)") ifun, 2.0
!old  end do
!old  do iact = 1, nact
!old     ifun = ncore + iact
!old     write(6, "(i5, e25.15)") ifun, dble(den1(iact, iact))
!old  end do
!old  do ifun = nocc + 1, nfun
!old     write(6, "(i5, e25.15)") ifun, 0.0
!old  end do  
!old
!old  !write(6, "('eigen-solutions of den1 (R):')")
!old  !do iact = 1, nact
!old  !   write(6, "(5x, i12)", advance = 'no') iact
!old  !end do
!old  !write(6, *)
!old  !do iact = 1, nact
!old  !   write(6, "(5x, f12.5)", advance = 'no') dble(den1(iact, iact))
!old  !end do
!old  !write(6, *)
!old  !do jact = 1, nact
!old  !   do iact = 1, nact
!old  !      write(6, "(i5, f12.5)", advance = 'no') jact, dble(umat(jact, iact))
!old  !   end do
!old  !   write(6, *)
!old  !end do
!old  !write(6, "('eigen-solutions of den1 (I):')")
!old  !do iact = 1, nact
!old  !   write(6, "(5x, i12)", advance = 'no') iact
!old  !end do
!old  !write(6, *)
!old  !do iact = 1, nact
!old  !   write(6, "(5x, f12.5)", advance = 'no') aimag(den1(iact, iact))
!old  !end do
!old  !write(6, *)
!old  !do jact = 1, nact
!old  !   do iact = 1, nact
!old  !      write(6, "(i5, f12.5)", advance = 'no') jact, aimag(umat(jact, iact))
!old  !   end do
!old  !   write(6, *)
!old  !end do
!old  ! DEBUG
!old
!old  nwfn(1:(nrad-1), 0:lmax1, 1:nfun) = czero
!old  do ifun = 1, ncore
!old     do l = abs(mval(ifun)), lmax1
!old        do irad = 1, nrad - 1
!old           nwfn(irad, l, ifun) = wfn(irad, l, ifun)
!old        end do
!old     end do
!old  end do
!old  do iact = 1, nact
!old     ifun = ncore + iact
!old     mi = mval(ifun)
!old
!old     do jact = 1, nact
!old        jfun = ncore + jact
!old        mj = mval(jfun)
!old
!old!!!     if (mi == mj) then
!old        m = 0
!old!!!     m = mi
!old        do l = abs(m), lmax1
!old           do irad = 1, nrad - 1
!old              nwfn(irad, l, ifun) = nwfn(irad, l, ifun) &
!old                                   + wfn(irad, l, jfun) * umat(jact, iact)
!old           end do
!old        end do
!old!!!     end if
!old     end do
!old  end do
!old  do ifun = nocc + 1, nfun
!old     do l = abs(mval(ifun)), lmax1
!old        do irad = 1, nrad - 1
!old           nwfn(irad, l, ifun) = wfn(irad, l, ifun)
!old        end do
!old     end do
!old  end do
!old
!old  deallocate(umat)

end subroutine hprod_norb
!######################################################################
!######################################################################
subroutine hprod_norb2(den1, umat, occn)

  use, intrinsic :: iso_c_binding
  use mod_const, only : czero
  use mod_rad, only : nrad
  use mod_sph, only : lmax1
  use mod_bas, only : mval
  use mod_ormas, only : ncore, nact, nocc, nfun

  implicit none
  complex(c_double_complex), intent(in) :: den1(1:nact, 1:nact)
  complex(c_double_complex), intent(out) :: umat(1:nact, 1:nact)
  real(c_double), intent(out) :: occn(1:nact)

  integer(c_long) :: iact, jact
  complex(c_double_complex), allocatable :: dtmp(:,:)

  if (nact > 0) then
     allocate(dtmp(1:nact, 1:nact))
     umat(1:nact, 1:nact) = czero
     dtmp(1:nact, 1:nact) = den1(1:nact, 1:nact)
     call lapack_zheev_dsc(nact, dtmp, umat)

     do iact = 1, nact
        occn(iact) = dble(dtmp(iact, iact))
     end do

     ! DEBUG
     write(6, "('hprod_norb2: 1RDM-Real')")
     write(6, "(10x)", advance = 'no')
     do jact = 1, nact
        write(6, "(i20)", advance = 'no') mval(ncore+jact)
     end do
     write(6, *)
     do iact = 1, nact
        write(6, "(2i5)", advance = 'no') iact, mval(ncore+iact)
        do jact = 1, nact
           write(6, "(e25.15)", advance = 'no') dble(den1(jact, iact))
        end do
        write(6, *)
     end do

     write(6, "('hprod_norb2: 1RDM-Imag')")
     write(6, "(10x)", advance = 'no')
     do jact = 1, nact
        write(6, "(i20)", advance = 'no') mval(ncore+jact)
     end do
     write(6, *)
     do iact = 1, nact
        write(6, "(2i5)", advance = 'no') iact, mval(ncore+iact)
        do jact = 1, nact
           write(6, "(e25.15)", advance = 'no') aimag(den1(jact, iact))
        end do
        write(6, *)
     end do

     write(6, "('hprod_norb2: natural occupation number')")
     write(6, "(30x)", advance = 'no')
     do jact = 1, nact
        write(6, "(i20)", advance = 'no') mval(ncore+jact)
     end do
     write(6, *)
     do iact = 1, nact
        write(6, "(2i5, e25.15)", advance = 'no') iact, mval(ncore+iact), occn(iact)
        do jact = 1, nact
           write(6, "(f10.5, f10.5)", advance = 'no') umat(jact, iact)
        end do
        write(6, *)
     end do

     deallocate(dtmp)
  end if

end subroutine hprod_norb2
!######################################################################
subroutine hprod_norb3(den1, umat, occn)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : mval
  use mod_rad, only : nrad
  use mod_sph, only : lmax1, mmax1
  use mod_const, only : zero, czero
  use mod_ormas, only : ncore, nact, nocc, nvir, nfun

  implicit none
  complex(c_double_complex), intent(inout) :: den1(1:nact, 1:nact)
  complex(c_double_complex), intent(out) :: umat(1:nact, 1:nact)
  real(c_double), intent(out) :: occn(1:nact)

  integer(c_long) :: iact, jact, ifun, jfun, irad, l, m, mi, mj, norb_msep, iorb, jorb, ntot, itot
  complex(c_double_complex) :: fac1, fac2
  integer(c_long), allocatable :: map_msep(:)
  complex(c_double_complex), allocatable :: den1_msep(:,:)
  complex(c_double_complex), allocatable :: umat_msep(:,:)

  if (nact <= 0) return
  occn(1:nact) = zero
  umat(1:nact, 1:nact) = czero

  allocate(map_msep(1:nact))

  ntot = 0
  do m = -mmax1, mmax1
     norb_msep = 0
     map_msep(1:nact) = 0
     do iact = 1, nact
        ifun = ncore + iact
        if (mval(ifun) == m) then
           norb_msep = norb_msep + 1
           map_msep(norb_msep) = iact
        end if
     end do

     if (norb_msep > 0) then
        allocate(den1_msep(1:norb_msep, 1:norb_msep))
        allocate(umat_msep(1:norb_msep, 1:norb_msep))

        do iorb = 1, norb_msep
           iact = map_msep(iorb)
           do jorb = 1, norb_msep
              jact = map_msep(jorb)
              den1_msep(jorb, iorb) = den1(jact, iact)
           end do
        end do
        !DEBUG
        write(6, "('hprod_norb3: m-separated 1RDM-Real: m = ', i5)") m
        write(6, "(10x)", advance = 'no')
        do jorb = 1, norb_msep
           write(6, "(i20)", advance = 'no') mval(ncore+map_msep(jorb))
        end do
        write(6, *)
        do iorb = 1, norb_msep
           write(6, "(2i5)", advance = 'no') iorb, mval(ncore+map_msep(iorb))
           do jorb = 1, norb_msep
              write(6, "(e25.15)", advance = 'no') dble(den1_msep(jorb, iorb))
           end do
           write(6, *)
        end do
        write(6, "('hprod_norb3: m-separated 1RDM-Imag: m = ', i5)") m
        write(6, "(10x)", advance = 'no')
        do jorb = 1, norb_msep
           write(6, "(i20)", advance = 'no') mval(ncore+map_msep(jorb))
        end do
        write(6, *)
        do iorb = 1, norb_msep
           write(6, "(2i5)", advance = 'no') iorb, mval(ncore+map_msep(iorb))
           do jorb = 1, norb_msep
              write(6, "(e25.15)", advance = 'no') aimag(den1_msep(jorb, iorb))
           end do
           write(6, *)
        end do
        !DEBUG

        umat_msep(1:norb_msep, 1:norb_msep) = czero
        call lapack_zheev_dsc(norb_msep, den1_msep, umat_msep)

        do iorb = 1, norb_msep
           iact = map_msep(iorb)
           occn(iact) = dble(den1_msep(iorb, iorb))
!!           itot = ntot + iorb
!!           occn(itot) = dble(den1_msep(iorb, iorb))
           do jorb = 1, norb_msep
              jact = map_msep(jorb)
              umat(jact, iact) = umat_msep(jorb, iorb)
!!              umat(jact, itot) = umat_msep(jorb, iorb)
           end do
        end do

        !DEBUG
        write(6, "('hprod_norb3: natural occupation number: m = ', i5)") m
        write(6, "(30x)", advance = 'no')
        do jorb = 1, norb_msep
           jact = map_msep(jorb)
           write(6, "(i10, i10)", advance = 'no') jact, mval(ncore+jact)
        end do
        write(6, *)
        do iorb = 1, norb_msep
           iact = map_msep(iorb)
           write(6, "(2i5, e25.15)", advance = 'no') iact, mval(ncore+iact), occn(iact)
           do jorb = 1, norb_msep
              jact = map_msep(jorb)
              write(6, "(f10.5, f10.5)", advance = 'no') umat(jact, iact)
           end do
           write(6, *)
        end do
        !DEBUG
        
        ntot = ntot + norb_msep
        deallocate(den1_msep)
        deallocate(umat_msep)
     end if
  end do

  ! DEBUG
  write(6, "('hprod_norb3: 1RDM-Real')")
  write(6, "(10x)", advance = 'no')
  do jact = 1, nact
     write(6, "(i20)", advance = 'no') mval(ncore+jact)
  end do
  write(6, *)
  do iact = 1, nact
     write(6, "(2i5)", advance = 'no') iact, mval(ncore+iact)
     do jact = 1, nact
        write(6, "(e25.15)", advance = 'no') dble(den1(jact, iact))
     end do
     write(6, *)
  end do

  write(6, "('hprod_norb3: 1RDM-Imag')")
  write(6, "(10x)", advance = 'no')
  do jact = 1, nact
     write(6, "(i20)", advance = 'no') mval(ncore+jact)
  end do
  write(6, *)
  do iact = 1, nact
     write(6, "(2i5)", advance = 'no') iact, mval(ncore+iact)
     do jact = 1, nact
        write(6, "(e25.15)", advance = 'no') aimag(den1(jact, iact))
     end do
     write(6, *)
  end do

  write(6, "('hprod_norb3: natural occupation number')")
  write(6, "(30x)", advance = 'no')
  do jact = 1, nact
     write(6, "(i10, i10)", advance = 'no') jact, mval(ncore+jact)
  end do
  write(6, *)
  do iact = 1, nact
     write(6, "(2i5, e25.15)", advance = 'no') iact, mval(ncore+iact), occn(iact)
     do jact = 1, nact
        write(6, "(f10.5, f10.5)", advance = 'no') umat(jact, iact)
     end do
     write(6, *)
  end do
  !DEBUG

  deallocate(map_msep)

end subroutine hprod_norb3
!######################################################################
subroutine hprod_norb4(den1, wfn, uwfn, occn)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : mval, lval
  use mod_rad, only : nrad
  use mod_sph, only : lmax1, mmax1
  use mod_const, only : zero, two, czero
  use mod_ormas, only : ncore, nact, nocc, nvir, nfun

  implicit none
  complex(c_double_complex), intent(inout) :: den1(1:nact, 1:nact)
  complex(c_double_complex), intent(in) :: wfn(1:(nrad-1), 0:lmax1, 1:*)
  complex(c_double_complex), intent(out) :: uwfn(1:(nrad-1), 0:lmax1, 1:*)
  real(c_double), intent(out) :: occn(1:*)

  integer(c_long) :: iact, jact, ifun, jfun, irad, l, m, mi, mj, norb_msep, iorb, jorb, ntot, itot
  complex(c_double_complex) :: fac1, fac2
  integer(c_long), allocatable :: map_msep(:)
  complex(c_double_complex), allocatable :: den1_msep(:,:)
  complex(c_double_complex), allocatable :: umat_msep(:,:)

  if (ncore > 0) then
     uwfn(1:(nrad-1), 0:lmax1, 1:ncore) = wfn(1:(nrad-1), 0:lmax1, 1:ncore)
  end if

  if (nact <= 0) return
  uwfn(1:(nrad-1), 0:lmax1, (ncore+1):nocc) = czero

  allocate(map_msep(1:nact))
  occn(1:nocc) = zero
  occn(1:ncore) = two

  ntot = 0
  do m = -mmax1, mmax1
     norb_msep = 0
     map_msep(1:nact) = 0
     do iact = 1, nact
        ifun = ncore + iact
        if (mval(ifun) == m) then
           norb_msep = norb_msep + 1
           map_msep(norb_msep) = iact
        end if
     end do
     if (norb_msep > 0) then
        allocate(den1_msep(1:norb_msep, 1:norb_msep))
        allocate(umat_msep(1:norb_msep, 1:norb_msep))

        do iorb = 1, norb_msep
           iact = map_msep(iorb)
           do jorb = 1, norb_msep
              jact = map_msep(jorb)
              den1_msep(jorb, iorb) = den1(jact, iact)
           end do
        end do
!        !DEBUG
!        write(6, "('hprod_norb4: m-separated 1RDM-Real: m = ', i5)") m
!        write(6, "(10x)", advance = 'no')
!        do jorb = 1, norb_msep
!           write(6, "(i20)", advance = 'no') mval(ncore+map_msep(jorb))
!        end do
!        write(6, *)
!        do iorb = 1, norb_msep
!           write(6, "(2i5)", advance = 'no') iorb, mval(ncore+map_msep(iorb))
!           do jorb = 1, norb_msep
!              write(6, "(e25.15)", advance = 'no') dble(den1_msep(jorb, iorb))
!           end do
!           write(6, *)
!        end do
!        write(6, "('hprod_norb4: m-separated 1RDM-Imag: m = ', i5)") m
!        write(6, "(10x)", advance = 'no')
!        do jorb = 1, norb_msep
!           write(6, "(i20)", advance = 'no') mval(ncore+map_msep(jorb))
!        end do
!        write(6, *)
!        do iorb = 1, norb_msep
!           write(6, "(2i5)", advance = 'no') iorb, mval(ncore+map_msep(iorb))
!           do jorb = 1, norb_msep
!              write(6, "(e25.15)", advance = 'no') aimag(den1_msep(jorb, iorb))
!           end do
!           write(6, *)
!        end do
!        !DEBUG
        umat_msep(1:norb_msep, 1:norb_msep) = czero
        call lapack_zheev_dsc(norb_msep, den1_msep, umat_msep)

        do iorb = 1, norb_msep
           iact = map_msep(iorb)
           ifun = ncore + iact
           occn(ncore+iact) = dble(den1_msep(iorb, iorb))
           do jorb = 1, norb_msep
              jact = map_msep(jorb)
              jfun = ncore + jact
              uwfn(1:(nrad-1), abs(m):lmax1, ifun) = &
              uwfn(1:(nrad-1), abs(m):lmax1, ifun) + &
               wfn(1:(nrad-1), abs(m):lmax1, jfun) * umat_msep(jorb, iorb)
           end do
        end do

        !DEBUG
        write(6, "('hprod_norb4: natural occupation number: m = ', i5)") m
        do iorb = 1, norb_msep
           iact = map_msep(iorb)
           write(6, "(2i5, e25.15)", advance = 'no') iact, mval(ncore+iact), occn(ncore+iact)
           write(6, *)
        end do
        !DEBUG
        
        ntot = ntot + norb_msep
        deallocate(den1_msep)
        deallocate(umat_msep)
     end if
  end do

  ! DEBUG
!  write(6, "('hprod_norb4: 1RDM-Real')")
!  write(6, "(10x)", advance = 'no')
!  do jact = 1, nact
!     write(6, "(i20)", advance = 'no') mval(ncore+jact)
!  end do
!  write(6, *)
!  do iact = 1, nact
!     write(6, "(2i5)", advance = 'no') iact, mval(ncore+iact)
!     do jact = 1, nact
!        write(6, "(e25.15)", advance = 'no') dble(den1(jact, iact))
!     end do
!     write(6, *)
!  end do
!
!  write(6, "('hprod_norb4: 1RDM-Imag')")
!  write(6, "(10x)", advance = 'no')
!  do jact = 1, nact
!     write(6, "(i20)", advance = 'no') mval(ncore+jact)
!  end do
!  write(6, *)
!  do iact = 1, nact
!     write(6, "(2i5)", advance = 'no') iact, mval(ncore+iact)
!     do jact = 1, nact
!        write(6, "(e25.15)", advance = 'no') aimag(den1(jact, iact))
!     end do
!     write(6, *)
!  end do

  !DEBUG
  write(6, "('hprod_norb4: natural occupation number')")
  do iact = 1, nact
     write(6, "(2i5, e25.15)") iact, mval(ncore+iact), occn(ncore+iact)
  end do
  !DEBUG

  deallocate(map_msep)

end subroutine hprod_norb4
!######################################################################
