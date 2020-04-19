!################################################################################
subroutine futil_gmatinv(n, thresh, a, inva)

  use, intrinsic :: iso_c_binding

  implicit none
  integer(c_long), intent(in) :: n
  real(c_double), intent(in) :: thresh
  complex(c_double_complex), intent(in) :: a(1:n, 1:n)
  complex(c_double_complex), intent(out) :: inva(1:n, 1:n)

  integer(c_long) :: i, j, k
  complex(c_double_complex), allocatable :: aa(:,:), raa(:,:)
  complex(c_double_complex), parameter :: czero = (0.d+0, 0.d+0)

  allocate(aa(1:n,1:n))
  allocate(raa(1:n,1:n))

  aa(1:n, 1:n) = czero
  do i = 1, n
     do j = 1, n
        do k = 1, n
           aa(j, i) = aa(j, i) + a(j, k) * conjg(a(i, k))
        end do
     end do
  end do
  call futil_matinv_reg(n, thresh, aa, raa)

  aa(1:n, 1:n) = a(1:n, 1:n)
  inva(1:n, 1:n) = czero
  do i = 1, n
     do j = 1, n
        do k = 1, n
           inva(j, i) = inva(j, i) + conjg(aa(k, j)) * raa(k, i)
        end do
     end do
  end do

  deallocate(raa)
  deallocate(aa)

end subroutine futil_gmatinv
!################################################################################
