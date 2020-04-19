!################################################################################
subroutine util_gmatinv(n, thresh, a, inva)
!
! input a and output inva are safely equivalenced.
!
  use, intrinsic :: iso_c_binding

  implicit none
  integer(c_int), intent(in) :: n
  real(c_double), intent(in) :: thresh
  complex(c_double_complex), intent(in) :: a(1:n, 1:n)
  complex(c_double_complex), intent(out) :: inva(1:n, 1:n)

  integer(c_int) :: i, j, k
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
  call util_matinv_reg(n, thresh, aa, raa)

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

end subroutine util_gmatinv
!################################################################################
subroutine util_gmatinv_real(n, thresh, a, inva)
!
! input a and output inva are safely equivalenced.
!
  use, intrinsic :: iso_c_binding

  implicit none
  integer(c_int), intent(in) :: n
  real(c_double), intent(in) :: thresh
  real(c_double), intent(in) :: a(1:n, 1:n)
  real(c_double), intent(out) :: inva(1:n, 1:n)

  integer(c_int) :: i, j, k
  real(c_double), allocatable :: aa(:,:), raa(:,:)

  allocate(aa(1:n,1:n))
  allocate(raa(1:n,1:n))

  aa(1:n, 1:n) = 0d0
  do i = 1, n
     do j = 1, n
        do k = 1, n
           aa(j, i) = aa(j, i) + a(j, k) * a(i, k)
        end do
     end do
  end do
  call util_matinv_reg_real(n, thresh, aa, raa)

  aa(1:n, 1:n) = a(1:n, 1:n)
  inva(1:n, 1:n) = 0d0
  do i = 1, n
     do j = 1, n
        do k = 1, n
           inva(j, i) = inva(j, i) + aa(k, j) * raa(k, i)
        end do
     end do
  end do

  deallocate(raa)
  deallocate(aa)

end subroutine util_gmatinv_real
!################################################################################
