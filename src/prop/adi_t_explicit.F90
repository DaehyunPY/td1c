!######################################################################
subroutine adi_t_explicit(tadi, wfn)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad, ndvr
  use mod_sph, only : lmax1
  use mod_bas, only : mval
  use mod_ormas, only : nfcore, nfun
  use mod_const, only : czero, runit

  implicit none
  complex(c_double_complex), intent(in) :: tadi(1:(2*ndvr+1), 1:(nrad-1), 0:lmax1)
  complex(c_double_complex), intent(inout) :: wfn(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), allocatable :: twfn(:,:,:)

  integer(c_int) :: ifun, l, dim, ld, lll, ull
  dim = nrad - 1
  ld = 2 * ndvr + 1
  allocate(twfn(1:(nrad-1), 0:lmax1, 1:nfun))

  !$omp parallel default(shared) private(lll, ull)
  !###########################
  call util_omp_disp(0, lmax1, lll, ull)
  do ifun = nfcore + 1, nfun
     do l = max(lll, abs(mval(ifun))), ull
        call ZGBMV('N', dim, dim, ndvr, ndvr, runit, tadi(1,1,l), &
             ld, wfn(1,l,ifun), 1, czero, twfn(1,l,ifun), 1)
        ! copy back
        wfn(1:(nrad-1), l, ifun) = twfn(1:(nrad-1), l, ifun)
     end do
  end do
  !###########################
  !$omp end parallel

  deallocate(twfn)

end subroutine adi_t_explicit
!######################################################################
