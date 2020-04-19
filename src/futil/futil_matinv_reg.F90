!################################################################################
subroutine futil_matinv_reg(n, thresh, a, inva)

  use, intrinsic :: iso_c_binding

  implicit none
  integer(c_long), intent(in) :: n
  real(c_double), intent(in) :: thresh
  complex(c_double_complex), intent(in) :: a(1:n, 1:n)
  complex(c_double_complex), intent(out) :: inva(1:n, 1:n)

  integer(c_long) :: ifun, jfun, kfun
  complex(c_double_complex) :: diag, invd
  complex(c_double_complex), allocatable :: uvec(:,:)
  complex(c_double_complex), allocatable :: work(:,:)
  complex(c_double_complex), parameter :: czero = (0.d+0, 0.d+0)
  complex(c_double_complex), parameter :: runit = (1.d+0, 0.d+0)

  allocate(uvec(1:n, 1:n))
  allocate(work(1:n, 1:n))

!td1d  call util_zcopy(n*n, a, 1, work, 1)
!td1d  call util_zcopy(n*n, czero, 0, inva, 1)
!td1d  call util_zcopy(n*n, czero, 0, uvec, 1)
  work(1:n, 1:n) = a(1:n, 1:n)
  inva(1:n, 1:n) = czero
  uvec(1:n, 1:n) = czero
  call futil_diag_comp(.true., n, work, uvec)
  
  do kfun = 1, n
     diag = work(kfun, kfun)
!     invd = runit / diag
     invd = diag / (diag * diag + thresh)
     do jfun = 1, n
        do ifun = 1, n
           inva(ifun, jfun) = inva(ifun, jfun) + uvec(ifun, kfun) * invd &
                                       & * conjg(uvec(jfun, kfun))
        end do
     end do
  end do

  deallocate(work)
  deallocate(uvec)

end subroutine futil_matinv_reg
!################################################################################
