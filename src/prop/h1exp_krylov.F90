!######################################################################
subroutine h1exp_krylov(isplit, igauge, lfield, dtime, maxcyc, thresh, &
     wfn, ntot, alph, beta, orb)
                        
  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad
  use mod_bas, only : nbas
  use mod_const, only : one, zero, czero

  implicit none
  integer(c_long), intent(in) :: isplit, igauge, maxcyc
  real(c_double), intent(in) :: lfield(3), dtime, thresh
  complex(c_double_complex), intent(in) :: wfn(1:nbas, 1:*)
  integer(c_long), intent(out) :: ntot
  real(c_double), intent(out) :: alph(1:maxcyc)
  real(c_double), intent(out) :: beta(1:maxcyc)
  complex(c_double_complex), intent(out) :: orb(1:nbas, 1:*)

  integer(c_long) :: icyc, jcyc
  real(c_double) :: oores  ! one over residue (1/res)
  real(c_double) :: error  ! error?
  complex(c_double_complex) :: tmp, zfield
  complex(c_double_complex), allocatable :: horb(:,:)
  real(c_double), external :: h1exp_krylov_error
!debug
!debug  real(c_double), parameter :: cutoff = 10.D+0
!debug

  allocate(horb(1:nbas, 1:1))

  ! ##### initialization #####
  zfield = lfield(3)
  alph(1:maxcyc) = zero
  beta(1:maxcyc) = zero
  orb(1:nbas, 1) = wfn(1:nbas, 1)

  ntot = 1

  do icyc = 1, maxcyc

     ! new sigma vectors
     horb(1:nbas, 1:1) = czero
     call hprod_tprod(orb(1, icyc), horb, 1, 1, 1, nrad - 1)
     if (isplit == 1 .and. igauge == 0) call hprod_zprod(zfield, orb(1, icyc), horb, 1, 1, 1, nrad - 1)
     if (isplit == 1 .and. igauge /= 0) call hprod_pzprod(zfield, orb(1, icyc), horb, 1, 1, 1, nrad - 1)

     ! new trial functions
     call zdotc_omp(nbas, orb(1, icyc), horb(1, 1), tmp)
     alph(icyc) = dble(tmp)
!debug     if (alph(icyc) > cutoff) then
!debug        write(6, "('h1exp_krylov: converged', i10, e20.10)") icyc - 1, error
!debug        ntot = icyc - 1
!debug        exit
!debug     end if

     if (icyc == 1) then
        horb(1:nbas, 1) = &
        horb(1:nbas, 1) - orb(1:nbas, icyc) * alph(icyc)
     else
        horb(1:nbas, 1) = &
        horb(1:nbas, 1) - orb(1:nbas, icyc)     * alph(icyc) &
                        - orb(1:nbas, icyc - 1) * beta(icyc - 1)
        ! ##### EXACT orthogonalization ######################
        do jcyc = icyc - 2, 1, -1
           call zdotc_omp(nbas, orb(1, jcyc), horb(1, 1), tmp)
           horb(1:nbas, 1) = horb(1:nbas, 1) - orb(1:nbas, jcyc) * tmp
        end do
        ! ####################################################
     end if

     call zdotc_omp(nbas, horb(1, 1), horb(1, 1), tmp)
     beta(icyc) = sqrt(dble(tmp))

     error = h1exp_krylov_error(icyc, alph, beta, dtime)
     write(6, "('h1exp_krylov: ', i10, 2f20.10, e20.10)") icyc, alph(icyc), beta(icyc), error

!    if (icyc == maxcyc .or. beta(icyc) < thresh) then
     if (error < thresh) then
        write(6, "('h1exp_krylov: converged', i10, e20.10)") icyc, error
        ntot = icyc
        exit
     else if (icyc == maxcyc) then
        write(6, "('h1exp_krylov: not converged', i10, e20.10)") icyc, error
        ntot = icyc
        exit
     else
        oores = one / beta(icyc)
        orb(1:nbas, icyc + 1) = horb(1:nbas, 1) * oores
!       orb(1:nbas, icyc + 1) = horb(1:nbas, 1)
     end if
  end do

  !DEBUG
  !do icyc = 1, ntot
  !   do jcyc = 1, icyc
  !      call zdotc_omp(nbas, orb(1, jcyc), orb(1, icyc), tmp)
  !      write(6, "(2i5, 2e20.10)") jcyc, icyc, tmp
  !   end do
  !end do
  !DEBUG

  deallocate(horb)

end subroutine h1exp_krylov
!######################################################################
subroutine h1exp_krylov1(isplit, igauge, lfield, maxcyc, maxvec, norb, mmap, wfn, &
     ntot, norm, orb, hmat)
                        
  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad
  use mod_bas, only : nbas
  use mod_ormas, only : nfun
  use mod_const, only : one, czero

  implicit none
  integer(c_long), intent(in) :: isplit, igauge
  real(c_double), intent(in) :: lfield(3)
  integer(c_long), intent(in) :: maxcyc, maxvec, norb, mmap(1:nfun)
  complex(c_double_complex), intent(in) :: wfn(1:nbas, 1:*)

  integer(c_long), intent(out) :: ntot
  real(c_double), intent(out) :: norm(1:maxvec)
  complex(c_double_complex), intent(out) :: orb(1:nbas, 1:*)
  complex(c_double_complex), intent(out) :: hmat(1:maxvec, 1:maxvec)

  real(c_double) :: rres, oores
  complex(c_double_complex) :: crres, zfield
  real(c_double), parameter :: thrsub = 1.D-10
  integer(c_long) :: ndone, nnew, nvec, icyc, ivec, jvec
  complex(c_double_complex), allocatable :: horb(:,:)

  zfield = lfield(3)
  allocate(horb(1:nbas, 1:norb))

  ! ##### initialization #####
  hmat(1:maxvec, 1:maxvec) = czero
  do ivec = 1, norb
     norm(ivec) = one
     orb(1:nbas, ivec) = wfn(1:nbas, mmap(ivec))
     write(6, "('h1exp_krylov: ', i10, f20.10)") ivec, norm(ivec)
  end do

  ntot = norb
  ndone = 0

  do icyc = 1, maxcyc

     ! new sigma vectors
     nvec = ntot - ndone
     horb(1:nbas, 1:nvec) = czero
     call hprod_tprod(orb(1, ndone+1), horb, 1, nvec, 1, nrad - 1)
     if (isplit == 1) then
        if (igauge == 0) then
           call hprod_zprod(zfield, orb(1, ndone+1), horb, 1, nvec, 1, nrad - 1)
        else
           call hprod_pzprod(zfield, orb(1, ndone+1), horb, 1, nvec, 1, nrad - 1)
        end if
     end if

     ! subspace hamiltonian, and
     ! new trial functions
     nnew = 0

     do ivec = 1, nvec
        do jvec = 1, ntot + nnew
           call zdotc_omp(nbas, orb(1, jvec), horb(1, ivec), hmat(jvec, ndone + ivec))
           hmat(ndone + ivec, jvec) = conjg(hmat(jvec, ndone + ivec))

           horb(1:nbas, ivec) = &
           horb(1:nbas, ivec) - orb(1:nbas, jvec) * hmat(jvec, ndone + ivec)
        end do

        call zdotc_omp(nbas, horb(1, ivec), horb(1, ivec), crres)
        rres = sqrt(abs(crres))

        if (icyc < maxcyc .and. rres > thrsub) then
           nnew = nnew + 1
           oores = one / rres
           norm(ntot + nnew) = rres
           orb(1:nbas, ntot + nnew) = horb(1:nbas, ivec) * oores
!          orb(1:nbas, ntot + nnew) = horb(1:nbas, ivec)
           write(6, "('h1exp_krylov: ', i10, f20.10)") ntot + nnew, norm(ntot + nnew)
        end if
     end do

     ndone = ntot
     ntot = ndone + nnew

  end do

  write(6, "('hmat-r:')")
  do ivec = 1, ntot
     do jvec = 1, ntot
        write(6, "(e14.5)", advance = 'no') dble(hmat(jvec, ivec))
     end do
     write(6, *)
  end do
  write(6, "('hmat-i:')")
  do ivec = 1, ntot
     do jvec = 1, ntot
        write(6, "(e14.5)", advance = 'no') aimag(hmat(jvec, ivec))
     end do
     write(6, *)
  end do

  deallocate(horb)

end subroutine h1exp_krylov1
!######################################################################
real(c_double) function h1exp_krylov_error(nsub, alph, beta, dtime)

  use, intrinsic :: iso_c_binding

  implicit none
  integer(c_long), intent(in) :: nsub
  real(c_double), intent(in) :: alph(1:*)
  real(c_double), intent(in) :: beta(1:*)
  real(c_double), intent(in) :: dtime

  integer(c_long) :: isub
  real(c_double) :: error

  error = (nsub - 1) * log(dtime)
  do isub = 1, nsub - 1
     error = error - log(dble(isub))
     error = error + log(beta(isub))
  end do
  error = error + log(beta(nsub))

  error = exp(error)
  error = error * error
  h1exp_krylov_error = error

end function h1exp_krylov_error
!######################################################################
