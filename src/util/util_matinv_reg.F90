!################################################################################
subroutine util_matinv_reg(n, thresh, a, inva)

  use, intrinsic :: iso_c_binding

  use mod_const, only : czero, runit

  implicit none
  integer(c_int), intent(in) :: n
  real(c_double), intent(in) :: thresh
  complex(c_double_complex), intent(in) :: a(1:n, 1:n)
  complex(c_double_complex), intent(out) :: inva(1:n, 1:n)

  integer(c_int) :: ifun, jfun, kfun
  complex(c_double_complex) :: diag, invd
  complex(c_double_complex), allocatable :: uvec(:,:)
  complex(c_double_complex), allocatable :: work(:,:)

  allocate(uvec(1:n, 1:n))
  allocate(work(1:n, 1:n))

  inva = 0d0
  work = 0d0
  uvec = 0d0
  call util_diag_comp(.true., n, work, uvec)
  
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

end subroutine util_matinv_reg
!################################################################################
!################################################################################
subroutine util_matinv_reg_real(n, thresh, a, inva)

  use, intrinsic :: iso_c_binding

  implicit none
  integer(c_int), intent(in) :: n
  real(c_double), intent(in) :: thresh
  real(c_double), intent(in) :: a(1:n, 1:n)
  real(c_double), intent(out) :: inva(1:n, 1:n)

  integer(c_int) :: ifun, jfun, kfun
  real(c_double) :: diag, invd
  real(c_double), allocatable :: uvec(:,:)
  real(c_double), allocatable :: work(:,:)

  allocate(uvec(1:n, 1:n))
  allocate(work(1:n, 1:n))

  work(1:n, 1:n) = a(1:n, 1:n)
  inva(1:n, 1:n) = 0.d+0
  uvec(1:n, 1:n) = 0.d+0
  call util_diag_real(.true., n, work, uvec)
  
  do kfun = 1, n
     diag = work(kfun, kfun)
     invd = diag / (diag * diag + thresh)
     do jfun = 1, n
        do ifun = 1, n
           inva(ifun, jfun) = inva(ifun, jfun) + uvec(ifun, kfun) * invd &
                                             & * uvec(jfun, kfun)
        end do
     end do
  end do

  deallocate(work)
  deallocate(uvec)

end subroutine util_matinv_reg_real
!################################################################################
