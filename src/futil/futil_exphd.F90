!################################################################################
subroutine futil_exphd(n, fac, hmat, exph)

  use, intrinsic :: iso_c_binding
  use mod_const, only : czero

  implicit none
  integer(c_long), intent(in) :: n
  complex(c_double_complex), intent(in) :: fac
  complex(c_double_complex), intent(in) :: hmat(1:n, 1:n)
  complex(c_double_complex), intent(out) :: exph(1:n, 1:n)

  integer(c_long) :: i, j, k
  complex(c_double_complex) :: expd
  complex(c_double_complex), allocatable :: htmp(:,:)
  complex(c_double_complex), allocatable :: uvec(:,:)

  allocate(htmp(1:n, 1:n))
  allocate(uvec(1:n, 1:n))

  htmp(1:n, 1:n) = hmat(1:n, 1:n)
  uvec(1:n, 1:n) = czero
  call futil_diag_comp(.false., n, htmp, uvec)

!DEBUG
!  write(6, "('eigenvalues of hmat:')")
!  do j = 1, n
!     write(6, "(f12.5)", advance = 'no') dble(htmp(j, j))
!  end do
!  write(6, *)
!DEBUG

  exph(1:n, 1:n) = czero
  do k = 1, n
     expd = exp(fac * htmp(k, k))
     do i = 1, n
        do j = 1, n
           exph(j, i) = exph(j, i) + uvec(j, k) * expd * conjg(uvec(i, k))
        end do
     end do
  end do

!debug  write(6, "('U_{1k}:')")
!debug  do k = 1, n
!debug     write(6, "(2e14.5)") uvec(1, k)
!debug  end do

  deallocate(uvec)
  deallocate(htmp)

end subroutine futil_exphd
!################################################################################
subroutine futil_exphd_cutoff(n, fac, cutoff, hmat, exph)

  use, intrinsic :: iso_c_binding
  use mod_const, only : czero

  implicit none
  integer(c_long), intent(in) :: n
  complex(c_double_complex), intent(in) :: fac
  real(c_double), intent(in) :: cutoff
  complex(c_double_complex), intent(in) :: hmat(1:n, 1:n)
  complex(c_double_complex), intent(out) :: exph(1:n, 1:n)

  logical :: valid
  integer(c_long) :: i, j, k
  complex(c_double_complex) :: expd
  complex(c_double_complex), allocatable :: htmp(:,:)
  complex(c_double_complex), allocatable :: uvec(:,:)

  allocate(htmp(1:n, 1:n))
  allocate(uvec(1:n, 1:n))

  htmp(1:n, 1:n) = hmat(1:n, 1:n)
  uvec(1:n, 1:n) = czero
  call futil_diag_comp(.false., n, htmp, uvec)

!DEBUG
  write(6, "('eigenvalues of hmat:')")
  do j = 1, n
     valid = dble(htmp(j, j)) <  cutoff
     write(6, "(i10, f12.5, l10)") j, dble(htmp(j, j)), valid
  end do
  write(6, *)
!DEBUG

  exph(1:n, 1:n) = czero
  do k = 1, n
     if (dble(htmp(j, j)) > cutoff) cycle

     expd = exp(fac * htmp(k, k))
     do i = 1, n
        do j = 1, n
           exph(j, i) = exph(j, i) + uvec(j, k) * expd * conjg(uvec(i, k))
        end do
     end do
  end do

!debug  write(6, "('U_{1k}:')")
!debug  do k = 1, n
!debug     write(6, "(2e14.5)") uvec(1, k)
!debug  end do

  deallocate(uvec)
  deallocate(htmp)

end subroutine futil_exphd_cutoff
!################################################################################
