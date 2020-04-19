!######################################################################
subroutine util_krylov(dtime, dim, maxcyc, thresh, ntot, norm1, alph, beta, vec, get_hvec)
                        
  use, intrinsic :: iso_c_binding
  use mod_const, only : one, zero, czero

  implicit none
  real(c_double), intent(in) :: dtime, thresh
  integer(c_long), intent(in) :: dim, maxcyc
  integer(c_long), intent(out) :: ntot
  real(c_double), intent(out) :: norm1
  real(c_double), intent(out) :: alph(1:maxcyc)
  real(c_double), intent(out) :: beta(1:maxcyc)
  complex(c_double_complex), intent(inout) :: vec(1:dim, 1:*)
  external get_hvec

  integer(c_long) :: icyc, jcyc
  real(c_double) :: oores  ! one over residue (1/res)
  real(c_double) :: error  ! error?
  complex(c_double_complex) :: tmp
  complex(c_double_complex), allocatable :: hvec(:)
  real(c_double), external :: util_krylov_error

  ! ##### initialization #####
  alph(1:maxcyc) = zero
  beta(1:maxcyc) = zero

  ! ##### normalize first vector #####
  call zdotc_omp(dim, vec(1,1), vec(1,1), tmp)
  norm1 = sqrt(abs(tmp))
  if (abs(norm1) < thresh) then
     write(6, "('util_krylov: norm1 = ', e20.10)") norm1
     ntot = 0
     return
  end if
  tmp = one / norm1
  vec(1:dim, 1) = vec(1:dim, 1) * tmp

  allocate(hvec(1:dim))

  ntot = 1

  do icyc = 1, maxcyc

     ! new sigma vectors
     hvec(1:dim) = czero
     call get_hvec(vec(1, icyc), hvec)

     ! new trial functions
     call zdotc_omp(dim, vec(1, icyc), hvec, tmp)
     alph(icyc) = dble(tmp)

     if (icyc == 1) then
        hvec(1:dim) = &
        hvec(1:dim) - vec(1:dim, icyc) * alph(icyc)
     else
        hvec(1:dim) = &
        hvec(1:dim) - vec(1:dim, icyc) * alph(icyc) &
                    - vec(1:dim, icyc - 1) * beta(icyc - 1)
        ! ##### EXACT orthogonalization ######################
        do jcyc = icyc - 2, 1, -1
           call zdotc_omp(dim, vec(1, jcyc), hvec, tmp)
           hvec(1:dim) = hvec(1:dim) - vec(1:dim, jcyc) * tmp
        end do
        ! ####################################################
     end if

     call zdotc_omp(dim, hvec, hvec, tmp)
     beta(icyc) = sqrt(dble(tmp))

     error = util_krylov_error(icyc, alph, beta, dtime)
     !DEBUG
     !write(6, "('util_krylov: ', i10, 2f20.10, e20.10)") icyc, alph(icyc), beta(icyc), error
     !DEBUG

!    if (icyc == maxcyc .or. beta(icyc) < thresh) then
     if (error < thresh) then
        write(6, "('util_krylov: converged', i10, e20.10)") icyc, error
        ntot = icyc
        exit
     else if (icyc == maxcyc) then
        write(6, "('util_krylov: not converged', i10, e20.10)") icyc, error
        ntot = icyc
        exit
     else
        oores = one / beta(icyc)
        vec(1:dim, icyc + 1) = hvec(1:dim) * oores
!       vec(1:dim, icyc + 1) = hvec(1:dim)
     end if
  end do

  !DEBUG
  !do icyc = 1, ntot
  !   do jcyc = 1, icyc
  !      call zdotc_omp(dim, vec(1, jcyc), vec(1, icyc), tmp)
  !      write(6, "(2i5, 2e20.10)") jcyc, icyc, tmp
  !   end do
  !end do
  !DEBUG

  deallocate(hvec)

end subroutine util_krylov
!######################################################################
real(c_double) function util_krylov_error(nsub, alph, beta, dtime)

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
  util_krylov_error = error

end function util_krylov_error
!######################################################################
