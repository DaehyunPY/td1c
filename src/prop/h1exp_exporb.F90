!################################################################################
subroutine h1exp_exporb(icomp, dtime, maxcyc, ntot, alph, beta, torb, wfn)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : nbas
  use mod_ormas, only : nfun
  use mod_const, only : zero, one, czero, runit, iunit

  implicit none
  integer(c_int), intent(in) :: icomp, maxcyc, ntot
  real(c_double), intent(in) :: dtime, alph(1:*), beta(1:*)
  complex(c_double_complex), intent(in) :: torb(1:nbas, 1:maxcyc)
  complex(c_double_complex), intent(out) :: wfn(1:nbas, 1:nfun)

  real(c_double), parameter :: cutoff = 10.D+0
  complex(c_double_complex) :: dtfac
  complex(c_double_complex), allocatable :: hsub(:,:)
  complex(c_double_complex), allocatable :: exph(:,:)
  integer(c_int) :: ivec, jvec, ibas, llb, ulb

  !DEBUG
  ! estimated error
  ! real(c_double) :: error
  ! integer(c_int), external :: util_ifact
  ! error = dtime ** (ntot - 1) / util_ifact(ntot - 1, ntot - 1)
  ! do ivec = 1, ntot
  !    error = error * beta(ivec)
  ! end do
  ! error = error ** 2
  ! write(6, "('estimated error = ', E20.10)") error
  !DEBUG

  if (icomp == 1) then
     dtfac = - dtime * iunit
  else
     dtfac = - dtime * runit
  end if
   
  allocate(hsub(1:ntot, 1:ntot))
  allocate(exph(1:ntot, 1:ntot))
  hsub(1:ntot, 1:ntot) = czero
  exph(1:ntot, 1:ntot) = czero
  do ivec = 1, ntot - 1
     hsub(ivec, ivec) = alph(ivec)
     hsub(ivec, ivec + 1) = beta(ivec)
     hsub(ivec + 1, ivec) = beta(ivec)
  end do
  hsub(ntot, ntot) = alph(ntot)

  !debug write(6, "('hsub:')")
  !debug do ivec = 1, ntot
  !debug    do jvec = 1, ntot
  !debug       write(6, "(e14.5)", advance = 'no') dble(hsub(jvec, ivec))
  !debug    end do
  !debug    write(6, *)
  !debug end do

! call futil_exphd_cutoff(ntot, dtfac, cutoff, hsub, exph)
  call futil_exphd(ntot, dtfac, hsub, exph)
! call futil_texphd(ntot, dtfac, hsub, exph)
! call futil_gexphd(ntot, dtfac, hsub, exph)

!debug  write(6, "('exph-r:')")
!debug  do ivec = 1, ntot
!debug     do jvec = 1, ntot
!debug        write(6, "(e14.5)", advance = 'no') dble(exph(jvec, ivec))
!debug     end do
!debug     write(6, *)
!debug  end do
!debug  write(6, "('exph-i:')")
!debug  do ivec = 1, ntot
!debug     do jvec = 1, ntot
!debug        write(6, "(e14.5)", advance = 'no') aimag(exph(jvec, ivec))
!debug     end do
!debug     write(6, *)
!debug  end do

  !$omp parallel default(shared) private(llb, ulb)
  call util_omp_disp(1, nbas, llb, ulb)

  wfn(llb:ulb, 1) = czero
  do jvec = 1, ntot
     do ibas = llb, ulb
        wfn(ibas, 1) =  wfn(ibas, 1) &
                     + torb(ibas, jvec) * exph(jvec, 1)
     end do
  end do
  !$omp end parallel

  deallocate(exph)
  deallocate(hsub)

end subroutine h1exp_exporb
!################################################################################
subroutine h1exp_exporb1(icomp, dtime, maxvec, ntot, norm, torb, hmat, norb, mmap, wfn)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : nbas
  use mod_ormas, only : nfun
  use mod_const, only : zero, one, czero, runit, iunit

  implicit none
  integer(c_int), intent(in) :: icomp, maxvec, ntot, norb, mmap(1:nfun)
  real(c_double), intent(in) :: dtime, norm(1:maxvec)
  complex(c_double_complex), intent(in) :: torb(1:nbas, 1:maxvec)
  complex(c_double_complex), intent(in) :: hmat(1:maxvec, 1:maxvec)
  complex(c_double_complex), intent(out) :: wfn(1:nbas, 1:nfun)

  complex(c_double_complex) :: dtfac
  complex(c_double_complex), allocatable :: hsub(:,:)
  complex(c_double_complex), allocatable :: exph(:,:)
  integer(c_int) :: ivec, jvec, ibas, iorb, ifun, llb, ulb

  if (icomp == 1) then
     dtfac = - dtime * iunit
  else
     dtfac = - dtime * runit
  end if
   
  allocate(hsub(1:ntot, 1:ntot))
  allocate(exph(1:ntot, 1:ntot))

  hsub(1:ntot, 1:ntot) = czero
  exph(1:ntot, 1:ntot) = czero
  do ivec = 1, ntot
     do jvec = 1, ntot
        hsub(jvec, ivec) = hmat(jvec, ivec)
!       hsub(jvec, ivec) = (hmat(jvec, ivec) / norm(jvec)) / norm(ivec)
     end do
  end do
  call futil_exphd(ntot, dtfac, hsub, exph)
! call futil_texphd(ntot, dtfac, hsub, exph)
! call futil_gexphd(ntot, dtfac, hsub, exph)

!  write(6, "('exph-r:')")
!  do ivec = 1, ntot
!     do jvec = 1, ntot
!        write(6, "(e14.5)", advance = 'no') dble(exph(jvec, ivec))
!     end do
!     write(6, *)
!  end do
!  write(6, "('exph-i:')")
!  do ivec = 1, ntot
!     do jvec = 1, ntot
!        write(6, "(e14.5)", advance = 'no') aimag(exph(jvec, ivec))
!     end do
!     write(6, *)
!  end do

  !$omp parallel default(shared) private(ifun, llb, ulb)
  call util_omp_disp(1, nbas, llb, ulb)
  do iorb = 1, norb
     ifun = mmap(iorb)
     wfn(llb:ulb, ifun) = czero
     do jvec = 1, ntot
        do ibas = llb, ulb
           wfn(ibas, ifun) =  wfn(ibas, ifun) &
                           + torb(ibas, jvec) * exph(jvec, iorb)
        end do
     end do
  end do
  !$omp end parallel

  deallocate(exph)
  deallocate(hsub)

end subroutine h1exp_exporb1
!######################################################################
