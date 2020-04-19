!################################################################################
subroutine futil_gexphd(n, fac, hmat, exph)

  use, intrinsic :: iso_c_binding

  implicit none
  integer(c_int), intent(in) :: n
  complex(c_double_complex), intent(in) :: fac
  complex(c_double_complex), intent(in) :: hmat(1:n, 1:n)
  complex(c_double_complex), intent(out) :: exph(1:n, 1:n)
  integer(c_int) :: i, j, k
  complex(c_double_complex) :: expd
  complex(c_double_complex), allocatable :: htmp(:,:)
  complex(c_double_complex), allocatable :: lvec(:,:), linv(:,:)
  complex(c_double_complex), allocatable :: rvec(:,:), rinv(:,:)
  real(c_double), parameter :: tiny = 1.D-20
  complex(c_double_complex), parameter :: czero = (0.d+0, 0.d+0)
  complex(c_double_complex), parameter :: runit = (1.d+0, 0.d+0)
  complex(c_double_complex), parameter :: iunit = (0.d+0, 1.d+0)

  allocate(htmp(1:n, 1:n))
  allocate(lvec(1:n, 1:n)); allocate(linv(1:n, 1:n))
  allocate(rvec(1:n, 1:n)); allocate(rinv(1:n, 1:n))

  htmp(1:n, 1:n) = hmat(1:n, 1:n)
  call futil_gdiag_comp(.false., n, htmp, lvec, rvec)
  call futil_gmatinv(n, tiny, lvec, linv)
  call futil_gmatinv(n, tiny, rvec, rinv)

  !DEBUG
  !DEBUG write(6, "('eigenvalues of hmat:')")
  !DEBUG do j = 1, n
  !DEBUG    write(6, "(i12, 2f12.8)") j, htmp(j, j)
  !DEBUG end do
  !DEBUG write(6, *)
  !DEBUG

  exph(1:n, 1:n) = czero
  do k = 1, n
     expd = exp(fac * htmp(k, k))
     !DEBUG
     !DEBUG write(6, "('futil_gexphd expd: ', i10, 2f20.10)") k, expd
     !DEBUG
     do i = 1, n
        do j = 1, n
           exph(j, i) = exph(j, i) + conjg(linv(k, j)) * expd * rinv(k, i)
        end do
     end do
  end do

  !DEBUG
  !DEBUG write(6, "('gexphd:')")
  !DEBUG do i = 1, n
  !DEBUG    do j = 1, n
  !DEBUG       write(6, "(2i10, 2f12.8)") i, j, exph(i, j)
  !DEBUG    end do
  !DEBUG end do
  !DEBUG write(6, *)
  !DEBUG

  deallocate(rvec); deallocate(rinv)
  deallocate(lvec); deallocate(linv)
  deallocate(htmp)

end subroutine futil_gexphd
!################################################################################
