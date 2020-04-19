!######################################################################
subroutine h1rat_prodci1(do_numer, do_denom, zfac, dtime, dimn, dimd, &
     ncoeff, dcoeff0, dcoeff1, h1type, eref, hdiag, int1e, int2e, cicin, cicout)

  use, intrinsic :: iso_c_binding
  use mod_const, only : zero, runit
  use mod_ormas, only : lcic, nact
  use mod_control, only : icomp, h1rat_maxcyc, h1rat_thresh

  implicit none
  logical(c_bool), intent(in) :: do_numer, do_denom
  complex(c_double_complex), intent(in) :: zfac
  real(c_double), intent(in) :: dtime, eref
  integer(c_long), intent(in) :: dimn, dimd, h1type
  complex(c_double_complex), intent(in) :: ncoeff(0:dimn)
  complex(c_double_complex), intent(in) :: dcoeff0(1:dimd), dcoeff1(1:dimd)
  complex(c_double_complex), intent(in) :: hdiag(1:lcic)
  complex(c_double_complex), intent(in) :: int1e(1:nact, 1:nact)
  complex(c_double_complex), intent(in) :: int2e(1:nact, 1:nact, 1:nact, 1:nact)
  complex(c_double_complex), intent(in) :: cicin(1:lcic)
  complex(c_double_complex), intent(inout) :: cicout(1:lcic)

  logical, parameter :: debug = .true.
  complex(c_double_complex), allocatable :: cic(:)
  complex(c_double_complex), allocatable :: dcic(:)
  complex(c_double_complex), allocatable :: tcic(:)
  complex(c_double_complex), allocatable :: d0inv(:)
  complex(c_double_complex) :: tfac, zfac1, cres
  integer(c_long) :: idim, ncyc, icyc, idet, ldet

  if (h1type .ne. 2) stop 'h1rat_prodci1: suport h1type = 2 only.'

  tfac = dtime
  ldet = lcic
  allocate(cic(1:ldet))
  allocate(dcic(1:ldet))
  allocate(tcic(1:ldet))
  allocate(d0inv(1:ldet))

!  call zcopy_omp(ldet, cicin, cic)
  call zcopy(ldet, cicin, 1, cic, 1)

  if (do_numer) then
!     call zcopy_omp(ldet, cic, dcic)
     call zcopy(ldet, cic, 1, dcic, 1)
     call zscal_omp(ldet, ncoeff(0), cic)
     do idim = 1, dimn
        call zclear_omp(ldet, tcic)
        call ormas_zhcicp(tfac, eref, int1e, int2e, dcic, tcic)
!        call zcopy_omp(ldet, tcic, dcic)
        call zcopy(ldet, tcic, 1, dcic, 1)
        call zaxpy_omp(ldet, ncoeff(idim), dcic, cic)
     end do
  end if

  if (do_denom) then
     do idim = 1, dimd
        if (abs(dcoeff1(idim)) < h1rat_thresh) cycle
     
        ! prepare 1/Q(0) = 1/(dcoeff0[i] + dcoeff1[i] * H0)
        !$omp parallel default(shared)
        !$omp do
        do idet = 1, lcic
           d0inv(idet) = runit / (dcoeff0(idim) + dcoeff1(idim) * hdiag(idet))
        end do
        !$omp end do
        !$omp end parallel

        ! 1/Q(0)
        !$omp parallel default(shared)
        !$omp do
        do idet = 1, lcic
           cic(idet) = d0inv(idet) * cic(idet)
        end do
        !$omp end do
        !$omp end parallel

        ! matrix iteration method
!        call zcopy_omp(ldet, cic, dcic)
        call zcopy(ldet, cic, 1, dcic, 1)
        do icyc = 1, h1rat_maxcyc
           ! convergence check
           ncyc = icyc
           call zdotc_omp(ldet, cic, cic, cres)
           if (sqrt(abs(cres)) < h1rat_thresh) exit
           if (debug) write(6, "('h1rat_prodci1: dim ', i5, ' cycle ', i5, ':', e20.5)") idim, icyc, sqrt(abs(cres))
     
           ! ######################################################
           ! - dcoeff1(i) * (H - H0)
           zfac1 = - dcoeff1(idim)
           call zclear_omp(ldet, tcic)
           call ormas_zhcicp(zfac, eref, int1e, int2e, cic, tcic)
           !$omp parallel default(shared)
           !$omp do
           do idet = 1, lcic
              tcic(idet) = tcic(idet) - d0inv(idet) * cic(idet)
           end do
           !$omp end do
           !$omp end parallel

           ! 1 / (dcoeff0[i] + dcoeff1[i] * H0)
           !$omp parallel default(shared)
           !$omp do
           do idet = 1, lcic
              cic(idet) = d0inv(idet) * tcic(idet)
              dcic(idet) = dcic(idet) + cic(idet)
           end do
           !$omp end do
           !$omp end parallel
           ! ######################################################
        end do
        write(6, "('h1rat_prodci1: dim ', i5, ' convd ', i5, ':', e20.5)") idim, ncyc, sqrt(abs(cres))

!        call zcopy_omp(ldet, dcic, cic)
        call zcopy(ldet, dcic, 1, cic, 1)
     end do
  end if

  !$omp parallel default(shared)
  !$omp do
  do idet = 1, lcic
     cicout(idet) = cicout(idet) + zfac * cic(idet)
  end do
  !$omp end do
  !$omp end parallel

  deallocate(d0inv)
  deallocate(tcic)
  deallocate(dcic)
  deallocate(cic)

end subroutine h1rat_prodci1
!######################################################################
