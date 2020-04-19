!######################################################################
subroutine adi_t_implicit(tadi, ipiv, wfn)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad, ndvr
  use mod_sph, only : lmax1
  use mod_bas, only : mval
  use mod_ormas, only : nfcore, nfun
  use mod_const, only : czero, runit

  implicit none
  complex(c_double_complex), intent(in) :: tadi(1:(3*ndvr+1), 1:(nrad-1), 0:lmax1)
  integer(c_long), intent(in) :: ipiv(1:(nrad-1), 0:lmax1)
  complex(c_double_complex), intent(inout) :: wfn(1:(nrad-1), 0:lmax1, 1:nfun)

  integer(c_long) :: ifun, l, dim, ld, info, lll, ull
  dim = nrad - 1
  ld = 3 * ndvr + 1

  !$omp parallel default(shared) private(lll, ull)
  !###########################
  call util_omp_disp(0, lmax1, lll, ull)
  do ifun = nfcore + 1, nfun
     do l = max(lll, abs(mval(ifun))), ull
        call ZGBTRS('N', dim, ndvr, ndvr, 1, tadi(1,1,l), ld, &
             ipiv(1,l), wfn(1,l,ifun), dim, info)
     end do
  end do
  !###########################
  !$omp end parallel

end subroutine adi_t_implicit
!######################################################################
