!######################################################################
subroutine h1rat_prodci2(do_numer, do_denom, zfac, dtime, dimn, dimd, &
     ncoeff, dcoeff0, dcoeff1, h1type, eref, hdiag, int1e, int2e, cicin, cicout)

  use, intrinsic :: iso_c_binding
  use mod_const, only : zero, czero, runit
  use mod_ormas, only : lcic, nact
  use mod_control, only : icomp, h1rat_maxcyc, h1rat_thresh
  use mod_h1rat, only: h1rat_eref, h1rat_int1e, h1rat_int2e

  implicit none
  logical(c_bool), intent(in) :: do_numer, do_denom
  complex(c_double_complex), intent(in) :: zfac
  real(c_double), intent(in) :: dtime, eref
  integer(c_int), intent(in) :: dimn, dimd, h1type
  complex(c_double_complex), intent(in) :: ncoeff(0:dimn)
  complex(c_double_complex), intent(in) :: dcoeff0(1:dimd), dcoeff1(1:dimd)
  complex(c_double_complex), intent(in) :: hdiag(1:lcic)
  complex(c_double_complex), intent(in) :: int1e(1:nact, 1:nact)
  complex(c_double_complex), intent(in) :: int2e(1:nact, 1:nact, 1:nact, 1:nact)
  complex(c_double_complex), intent(in) :: cicin(1:lcic)
  complex(c_double_complex), intent(inout) :: cicout(1:lcic)

  real(c_double) :: krylov_thresh
  real(c_double) :: norm1
  real(c_double) :: alph(1:h1rat_maxcyc)
  real(c_double) :: beta(1:h1rat_maxcyc)
  real(c_double), allocatable :: hmat(:,:)
  real(c_double), allocatable :: uvec(:,:)
  complex(c_double_complex) :: val, tmp
  complex(c_double_complex), allocatable :: func(:)
  complex(c_double_complex), allocatable :: cvec(:,:)
  integer(c_int) :: idim, nkry, isub, jsub, idet, ldet
  logical, parameter :: debug = .false.
  external :: h1rat_hcic1, h1rat_hcic

! krylov_thresh = 1.0D-15
  krylov_thresh = h1rat_thresh

  if (.not. do_denom) then
     call h1rat_prodci(do_numer, do_denom, zfac, dtime, dimn, dimd, &
          ncoeff, dcoeff0, dcoeff1, h1type, eref, hdiag, int1e, int2e, cicin, cicout)
     return
  end if

  ldet = lcic
  allocate(h1rat_int1e(1:nact, 1:nact))
  allocate(h1rat_int2e(1:nact, 1:nact, 1:nact, 1:nact))
  allocate(cvec(1:lcic, 1:h1rat_maxcyc))

!  call zcopy_omp(ldet, cicin, cvec(1,1))
  call zcopy(ldet, cicin, 1, cvec(1,1), 1)
  h1rat_eref = eref
  h1rat_int1e(1:nact, 1:nact) = int1e(1:nact, 1:nact)
  h1rat_int2e(1:nact, 1:nact, 1:nact, 1:nact) = int2e(1:nact, 1:nact, 1:nact, 1:nact)
  if (h1type == 1) then
     call util_krylov(dtime, ldet, h1rat_maxcyc, krylov_thresh, nkry, norm1, alph, beta, cvec, h1rat_hcic1)
  else 
     call util_krylov(dtime, ldet, h1rat_maxcyc, krylov_thresh, nkry, norm1, alph, beta, cvec, h1rat_hcic)
  end if

  if (nkry > 0) then

     allocate(hmat(1:nkry, 1:nkry))
     allocate(uvec(1:nkry, 1:nkry))
     allocate(func(1:nkry))
     hmat(1:nkry, 1:nkry) = zero
     uvec(1:nkry, 1:nkry) = zero
     do isub = 1, nkry
        hmat(isub, isub) = alph(isub)
        if (isub < nkry) hmat(isub, isub + 1) = beta(isub)
        if (isub >    1) hmat(isub, isub - 1) = beta(isub-1)
     end do
     if (debug) then
        write(6, "('Lanczos vectors')")
        do isub = 1, h1rat_maxcyc
           write(6, "(i5, 2f15.8)") isub, alph(isub), beta(isub)
        end do
        write(6, "('CI Hamiltonian in ', i5, ' dimensional Krylov space:')") nkry
        do isub = 1, nkry
           do jsub = 1, nkry
              write(6, "(f15.8)", advance = 'no') hmat(jsub, isub)
           end do
           write(6, *)
        end do
     end if

     call lapack_dsyev(nkry, hmat, uvec)

     if (debug) write(6, "('Eigenvalues and operator diagonals in the Krylov space:')")

     func(1:nkry) = czero
     do isub = 1, nkry
        val = runit
        ! numerator
        if (do_numer) then
           tmp = runit
           val = ncoeff(0) * val
           do idim = 1, dimn
              tmp = dtime * hmat(isub, isub) * tmp
              val = val + ncoeff(idim) * tmp
           end do
        end if
        ! denominator
        if (do_denom) then
           do idim = 1, dimd
              if (abs(dcoeff1(idim)) < h1rat_thresh) cycle
              val = val / (dcoeff0(idim) + dcoeff1(idim) * hmat(isub, isub))
           end do
        end if
        if (debug) write(6, "(i5, 3f15.8)") isub, hmat(isub, isub), val
        do jsub = 1, nkry
           func(jsub) = func(jsub) + uvec(jsub, isub) * val * uvec(1, isub) * norm1
        end do
     end do

     do isub = 1, nkry
        !$omp parallel default(shared)
        !$omp do
        do idet = 1, lcic
           cicout(idet) = cicout(idet) + zfac * cvec(idet, isub) * func(isub)
        end do
        !$omp end do
        !$omp end parallel
     end do

     deallocate(func)
     deallocate(uvec)
     deallocate(hmat)

  end if

  deallocate(cvec)
  deallocate(h1rat_int2e)
  deallocate(h1rat_int1e)

end subroutine h1rat_prodci2
!######################################################################
subroutine h1rat_hcic(cic, hcic)

  use, intrinsic :: iso_c_binding
  use mod_const, only : zero, czero
  use mod_ormas, only : lcic
  use mod_h1rat, only: h1rat_eref, h1rat_int1e, h1rat_int2e

  implicit none
  complex(c_double_complex), intent(in) :: cic(1:*)
  complex(c_double_complex), intent(out) :: hcic(1:*)

  hcic(1:lcic) = czero
  call ormas_hcic(h1rat_int1e, h1rat_int2e, cic, hcic, h1rat_eref)

end subroutine h1rat_hcic
!######################################################################
subroutine h1rat_hcic1(cic, hcic)

  use, intrinsic :: iso_c_binding
  use mod_const, only : czero
  use mod_h1rat, only: h1rat_int1e
  use mod_ormas, only : lcic

  implicit none
  complex(c_double_complex), intent(in) :: cic(1:*)
  complex(c_double_complex), intent(out) :: hcic(1:*)

  hcic(1:lcic) = czero
  call ormas_hcic1(h1rat_int1e, cic, hcic)

end subroutine h1rat_hcic1
!######################################################################
