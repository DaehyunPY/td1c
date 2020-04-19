!################################################################################
subroutine lapack_zhesv(n, nrhs, amat, bvec)

  use, intrinsic :: iso_c_binding

  implicit none
  integer(c_int), intent(in) :: n, nrhs
  complex(c_double_complex), intent(inout) :: amat(1:n, 1:n)
  complex(c_double_complex), intent(out) :: bvec(1:n, 1:nrhs)

  integer(c_int) :: info, lwork, i
  complex(c_double_complex) :: clwork
  integer(c_int), allocatable :: ipiv(:)
  complex(c_double_complex), allocatable :: work(:)

  allocate(ipiv(1:n))
  call zhesv("u", n, nrhs, amat, n, ipiv, bvec, n, clwork, -1, info)
  lwork = int(clwork)
  allocate(work(lwork))
  call zhesv("u", n, nrhs, amat, n, ipiv, bvec, n, work, lwork, info)

  if (info /= 0) then
     write(6, "('zhesv: info = ', i20)") info
     stop 'error in lapack_zhesv.'
  end if

  deallocate(work)
  deallocate(ipiv)

end subroutine lapack_zhesv
!################################################################################
subroutine lapack_zgesv(n, nrhs, amat, bvec)

  use, intrinsic :: iso_c_binding

  implicit none
  integer(c_int), intent(in) :: n, nrhs
  complex(c_double_complex), intent(inout) :: amat(1:n, 1:n)
  complex(c_double_complex), intent(out) :: bvec(1:n, 1:nrhs)

  integer(c_int) :: info, i
  integer(c_int), allocatable :: ipiv(:)

  allocate(ipiv(1:n))
  call zgesv(n, nrhs, amat, n, ipiv, bvec, n, info)
! call ZGESV(N, NRHS, A,  LDA, IPIV,  B, LDB, INFO)

  if (info /= 0) then
     write(6, "('zhesv: info = ', i20)") info
     stop 'error in lapack_zhesv.'
  end if

  deallocate(ipiv)

end subroutine lapack_zgesv
!################################################################################
