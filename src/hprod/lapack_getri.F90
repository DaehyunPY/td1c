!################################################################################
subroutine lapack_zgetri(n, a, u)

  use, intrinsic :: iso_c_binding

  implicit none
  integer(c_int), intent(in) :: n
  complex(c_double_complex), intent(inout) :: a(1:n, 1:n)
  complex(c_double_complex), intent(out) :: u(1:n, 1:n)

  integer(c_int) :: nf90, info, lwork, i, len
  complex(c_double_complex) :: clwork

  logical, parameter :: dsc = .false.
  complex(c_double_complex), allocatable :: atmp(:, :)
  complex(c_double_complex), allocatable :: work(:)
  real(c_double), allocatable :: eig(:)
  real(c_double), allocatable :: rwork(:)
!debug
!debug  integer(c_int) :: j, k
!debug

  nf90 = n
  len = max(nf90*nf90, 3*nf90-2)
  allocate(rwork(len))
  allocate(atmp(nf90, nf90))
  allocate(eig(nf90))
  atmp(1:nf90, 1:nf90) = a(1:nf90, 1:nf90)

!debug
!debugwrite(6, "('A: A')")
!debugwrite(6, "(2f20.10)") atmp(1:nf90, 1:nf90)
!debug
  call zgetri("v", "u", nf90, atmp, nf90, eig, clwork, -1, rwork, info)
  lwork = int(clwork)

  allocate(work(lwork))
  call zgetri("v", "u", nf90, atmp, nf90, eig, work, lwork, rwork, info)
!debug
!debugwrite(6, "('A: A = U * a * U+')")
!debugdo i = 1, nf90
!debug   do j = 1, nf90
!debug      a(i, j) = (0.d+0, 0.d+0)
!debug      do k = 1, nf90
!debug         a(i, j) = a(i, j) + atmp(i, k) * eig(k) * conjg(atmp(j, k))
!debug      end do
!debug   end do
!debugend do
!debugwrite(6, "(2f20.10)") a(1:nf90, 1:nf90)
!debug

  if (info /= 0) then
     write(6, "('util_diag: info = ', i20)") info
     stop 'error in util_diag_comp.'
  else
     u(1:nf90, 1:nf90) = (0.d+0, 0.d+0)
     a(1:nf90, 1:nf90) = (0.d+0, 0.d+0)
     if (dsc) then
        do i = 1, nf90
           a(i, i) = eig(nf90 - i + 1)
           u(1:nf90, i) = atmp(1:nf90, nf90 - i + 1)
        end do
     else
        do i = 1, nf90
           a(i, i) = eig(i)
           u(1:nf90, i) = atmp(1:nf90, i)
        end do
     end if
  end if

!debug
!write(6,"('diag: eigenvalues:')")
!write(6,"(F20.10)") eig(1:nf90)
!debug

  deallocate(work)
  deallocate(eig)
  deallocate(atmp)
  deallocate(rwork)

end subroutine lapack_zgetri
!################################################################################
