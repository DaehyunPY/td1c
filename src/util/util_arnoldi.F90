!######################################################################
subroutine util_arnoldi(dtime, dim, maxcyc, thresh, ntot, hmat, vec, get_hvec)
                        
  use, intrinsic :: iso_c_binding
  use mod_const, only : one, zero, czero, runit

  implicit none
  real(c_double), intent(in) :: dtime, thresh
  integer(c_int), intent(in) :: dim, maxcyc
  integer(c_int), intent(out) :: ntot
  complex(c_double_complex), intent(out) :: hmat(1:maxcyc, 1:maxcyc)
  complex(c_double_complex), intent(inout) :: vec(1:dim, 1:*)
  external get_hvec

  integer(c_int) :: icyc, jcyc
  real(c_double) :: error, beta, rbeta
  complex(c_double_complex) :: tmp
  complex(c_double_complex), allocatable :: hvec(:)
  real(c_double), external :: util_arnoldi_error

  allocate(hvec(1:dim))

  ! ##### initialization #####
  hmat(1:maxcyc, 1:maxcyc) = czero

  ! ##### normalize first vector #####
  call zdotc_omp(dim, vec(1,1), vec(1,1), tmp)
  tmp = runit / sqrt(abs(tmp))
  vec(1:dim, 1) = vec(1:dim, 1) * tmp

  ntot = 1

  do icyc = 1, maxcyc

     ! new sigma vectors
     hvec(1:dim) = czero
     call get_hvec(vec(1,icyc), hvec)

     ! orthogonalization
     do jcyc = 1, icyc
        call zdotc_omp(dim, vec(1,jcyc), hvec, tmp)
        hmat(jcyc, icyc) = tmp
        hvec(1:dim) = hvec(1:dim) - vec(1:dim, jcyc) * tmp
     end do
        
     ! residual
     call zdotc_omp(dim, hvec, hvec, tmp)
     beta = sqrt(abs(tmp))

     error = util_arnoldi_error(maxcyc, icyc, hmat, beta, dtime)
     write(6, "('util_arnoldi: ', i10, e20.10)") icyc, error

     if (error < thresh) then
        write(6, "('util_arnoldi: converged', i10, e20.10)") icyc, error
        ntot = icyc
        exit
     else if (icyc == maxcyc) then
        write(6, "('util_arnoldi: not convd', i10, e20.10)") icyc, error
        ntot = icyc
        exit
     else
        rbeta = one / beta
        hmat(icyc+1, icyc) = beta
        vec(1:dim, icyc+1) = hvec(1:dim) * rbeta
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

end subroutine util_arnoldi
!######################################################################
real(c_double) function util_arnoldi_error(maxcyc, nsub, hmat, beta, dtime)

  use, intrinsic :: iso_c_binding

  implicit none
  integer(c_int), intent(in) :: maxcyc, nsub
  complex(c_double_complex), intent(in) :: hmat(1:maxcyc, 1:maxcyc)
  real(c_double), intent(in) :: beta, dtime

  integer(c_int) :: isub
  real(c_double) :: error

  error = (nsub - 1) * log(dtime)
  do isub = 1, nsub - 1
     error = error - log(dble(isub))
     error = error + log(dble(hmat(isub+1, isub)))
  end do
  error = error + log(beta)

  error = exp(error)
  error = error * error
  util_arnoldi_error = error

end function util_arnoldi_error
!######################################################################
