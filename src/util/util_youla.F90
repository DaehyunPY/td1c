!################################################################################
subroutine util_youla_ceig(dim, amat, aeig, uvec)
!
  use, intrinsic :: iso_c_binding
!
! complex skew-symmetric matrix is factorized to the 
! direct sum of 2x2 real skew-symmetric factors.
! input  n: dimension
!        a: skew-symmetric matrix
! output a: diaga
!           diaga(2i-1,2i) =  sqrt(diag(aa+))
!           diaga(2i,2i-1) = -sqrt(diag(aa+))
!        u: unitary transformation matrix
!           diaga = u^T * a * u
!           a = v * diaga * v^T, v = conjg(u)
!
  implicit none
  integer(c_int), intent(in) :: dim
  complex(c_double_complex), intent(in) :: amat(1:dim, 1:dim)
  complex(c_double_complex), intent(out) :: aeig(1:*)
  complex(c_double_complex), intent(out) :: uvec(1:dim, 1:dim)

  integer(c_int) :: i, j, k, l, ii, dim2
  complex(c_double_complex) :: tmp
  complex(c_double_complex), allocatable :: atmp(:,:)
  complex(c_double_complex), allocatable :: vvec(:,:)
  integer(c_int), parameter :: iprint = 0

  allocate(atmp(1:dim, 1:dim))
  allocate(vvec(1:dim, 1:dim))

  vvec(1:dim, 1:dim) = (0.d+0, 0.d+0)
  atmp(1:dim, 1:dim) = (0.d+0, 0.d+0)
  do i = 1, dim
     do j = 1, dim
        do k = 1, dim
           atmp(i, j) = atmp(i, j) + amat(i, k) * conjg(amat(j, k))
        end do
     end do
  end do
  call util_diag_comp(.true., dim, atmp, vvec)

  do i = 1, dim
     do j = 1, dim
        uvec(i, j) = conjg(vvec(i, j))
     end do
  end do

  ! adiag = U^T a U
  atmp(1:dim, 1:dim) = (0.d+0, 0.d+0)
  do i = 1, dim
     do j = 1, dim
        do k = 1, dim
           do l = 1, dim
              atmp(i, j) = atmp(i, j) + uvec(k, i) * amat(k, l) * uvec(l, j)
           end do
        end do
     end do
  end do

  if (iprint > 4) then
     write(6, "('Youla decomposition of a. real-part:')")
     do i = 1, dim
        do j = 1, dim
           write(6, "(f12.5)", advance = 'no') dble(atmp(i, j))
        end do
        write(6, *)
     end do
     write(6, "('Youla decomposition of a. imag-part:')")
     do i = 1, dim
        do j = 1, dim
           write(6, "(f12.5)", advance = 'no') aimag(atmp(i, j))
        end do
        write(6, *)
     end do
  end if

  dim2 = dim / 2
  do ii = 1, dim2
     i = 2 * ii - 1
     aeig(ii) = atmp(i, i+1)
  end do

  write(6, "('util_youla_ceig. eigenvalues:')")
  do ii = 1, dim2
     write(6, "(i10, 2f16.8, e16.4)") ii, aeig(ii), abs(aeig(ii))
  end do
!debug  write(6, "('util_youla. eigenvectors:')")
!debug  do ii = 1, dim2
!debug     i = 2 * ii - 1
!debug     write(6, "(i16)", advance = 'no') i
!debug  end do
!debug  write(6, *)
!debug  do j = 1, dim
!debug     do ii = 1, dim2
!debug        i = 2 * ii - 1
!debug        write(6, "(e16.4)", advance = 'no') dble(uvec(j, i))
!debug     end do
!debug     write(6, *)
!debug  end do
!debug  do ii = 1, dim2
!debug     i = 2 * ii
!debug     write(6, "(i16)", advance = 'no') i
!debug  end do
!debug  write(6, *)
!debug  do j = 1, dim
!debug     do ii = 1, dim2
!debug        i = 2 * ii
!debug        write(6, "(e16.4)", advance = 'no') dble(uvec(j, i))
!debug     end do
!debug     write(6, *)
!debug  end do

  deallocate(vvec)
  deallocate(atmp)

!debug
!  stop "for debug in util_youla"
!debug

end subroutine util_youla_ceig
!################################################################################
!################################################################################
subroutine util_youla_reig(dim, amat, aeig, uvec)
!
  use, intrinsic :: iso_c_binding
!
! complex skew-symmetric matrix is factorized to the 
! direct sum of 2x2 real skew-symmetric factors.
! input  n: dimension
!        a: skew-symmetric matrix
! output a: diaga
!           diaga(2i-1,2i) =  sqrt(diag(aa+))
!           diaga(2i,2i-1) = -sqrt(diag(aa+))
!        u: unitary transformation matrix
!           diaga = u^T * a * u
!           a = v * diaga * v^T, v = conjg(u)
!
  implicit none
  integer(c_int), intent(in) :: dim
  complex(c_double_complex), intent(in) :: amat(1:dim, 1:dim)
  complex(c_double_complex), intent(out) :: aeig(1:*)
  complex(c_double_complex), intent(out) :: uvec(1:dim, 1:dim)

  integer(c_int) :: i, j, k, l, ii, dim2
  complex(c_double_complex) :: tmp
  complex(c_double_complex), allocatable :: atmp(:,:)
  complex(c_double_complex), allocatable :: vvec(:,:)

  allocate(atmp(1:dim, 1:dim))
  allocate(vvec(1:dim, 1:dim))

  vvec(1:dim, 1:dim) = (0.d+0, 0.d+0)
  atmp(1:dim, 1:dim) = (0.d+0, 0.d+0)
  do i = 1, dim
     do j = 1, dim
        do k = 1, dim
           atmp(i, j) = atmp(i, j) + amat(i, k) * conjg(amat(j, k))
        end do
     end do
  end do
  call util_diag_comp(.true., dim, atmp, vvec)

  do i = 1, dim
     do j = 1, dim
        uvec(i, j) = conjg(vvec(i, j))
     end do
  end do

  ! adiag = U^T a U
  atmp(1:dim, 1:dim) = (0.d+0, 0.d+0)
  do i = 1, dim
     do j = 1, dim
        do k = 1, dim
           do l = 1, dim
              atmp(i, j) = atmp(i, j) + uvec(k, i) * amat(k, l) * uvec(l, j)
           end do
        end do
     end do
  end do

  dim2 = dim / 2
  do ii = 1, dim2
     i = 2 * ii - 1
     tmp = abs(atmp(i, i+1))
     tmp = atmp(i, i+1) / tmp
     tmp = (1.d+0, 0.d+0) / sqrt(tmp)
     do j = 1, dim
        uvec(j, i)   = uvec(j, i)   * tmp
        uvec(j, i+1) = uvec(j, i+1) * tmp
     end do
  end do

!debug
!  vvec(1:dim, 1:dim) = (0.d+0, 0.d+0)
!  do i = 1, dim
!     do j = 1, dim
!        do k = 1, dim
!           do l = 1, dim
!              vvec(i, j) = vvec(i, j) + uvec(k, i) * amat(k, l) * uvec(l, j)
!           end do
!        end do
!     end do
!  end do
!  write(6, "('Youla decomposition of a. real-part:')")
!  do i = 1, dim
!     do j = 1, dim
!        write(6, "(f12.5)", advance = 'no') dble(vvec(i, j))
!     end do
!     write(6, *)
!  end do
!  write(6, "('Youla decomposition of a. imag-part:')")
!  do i = 1, dim
!     do j = 1, dim
!        write(6, "(f12.5)", advance = 'no') aimag(vvec(i, j))
!     end do
!     write(6, *)
!  end do
!debug

  do ii = 1, dim2
     i = 2 * ii - 1
     aeig(ii) = abs(atmp(i, i+1))
  end do

  write(6, "('util_youla_reig. eigenvalues:')")
  do ii = 1, dim2
     write(6, "(i10, 2f16.8, e16.4)") ii, aeig(ii), abs(aeig(ii))
  end do

  deallocate(vvec)
  deallocate(atmp)

end subroutine util_youla_reig
!################################################################################
!################################################################################
subroutine util_youla_real(dim, amat, aeig, uvec, vvec)
!
  use, intrinsic :: iso_c_binding
!
! real skew-symmetric matrix is factorized into the 
! direct sum of 2x2 real skew-symmetric factors.
! input   dim: dimension
!        amat: skew-symmetric matrix
! output aeig: diaga
!              diaga(2i-1,2i) =  sqrt(diag(aa+))
!              diaga(2i,2i-1) = -sqrt(diag(aa+))
!        uvec: unitary transformation matrix
!        vvec: unitary transformation matrix
!              diaga = u^T * a * v
!              a = u * diaga * v^T
!
  implicit none
  integer(c_int), intent(in) :: dim
  real(c_double), intent(in) :: amat(1:dim, 1:dim)
  real(c_double), intent(out) :: aeig(1:*)
  real(c_double), intent(out) :: uvec(1:dim, 1:dim)
  real(c_double), intent(out) :: vvec(1:dim, 1:dim)

  integer(c_int) :: i, j, k, l, ii, dim2
  real(c_double) :: tmp
  real(c_double), allocatable :: atmp(:,:)

  allocate(atmp(1:dim, 1:dim))

  vvec(1:dim, 1:dim) = 0.d+0
  atmp(1:dim, 1:dim) = 0.d+0
  do i = 1, dim
     do j = 1, dim
        do k = 1, dim
           atmp(i, j) = atmp(i, j) + amat(i, k) * amat(j, k)
        end do
     end do
  end do
  call util_diag_real(.true., dim, atmp, vvec)

  do i = 1, dim
     do j = 1, dim
        uvec(i, j) = vvec(i, j)
     end do
  end do

!debug  write(6, "('util_youla_real: A')")
!debug  do i = 1, dim
!debug     do j = 1, dim
!debug        write(6, "(f20.10)", advance = 'no') amat(i, j)
!debug     end do
!debug     write(6, *)
!debug  end do
!debug
!debug  write(6, "('util_youla_real. eigenvalues of aa^t:')")
!debug  do i = 1, dim
!debug     write(6, "(i10, f20.10)") i, atmp(i, i)
!debug  end do

  ! adiag = U^T a U
  atmp(1:dim, 1:dim) = 0.d+0
  do i = 1, dim
     do j = 1, dim
        do k = 1, dim
           do l = 1, dim
              atmp(i, j) = atmp(i, j) + uvec(k, i) * amat(k, l) * uvec(l, j)
           end do
        end do
     end do
  end do

!debug  write(6, "('util_youla_real: a = Vt * A * V')")
!debug  do i = 1, dim
!debug     do j = 1, dim
!debug        write(6, "(f20.10)", advance = 'no') atmp(i, j)
!debug     end do
!debug     write(6, *)
!debug  end do

  dim2 = dim / 2
  do ii = 1, dim2
     i = 2 * ii - 1
     if (atmp(i, i+1) < 0.d+0) then
        do j = 1, dim
           uvec(j, i)   = -uvec(j, i)
           uvec(j, i+1) = -uvec(j, i+1)
        end do
     end if
  end do

  ! adiag = U^T a V
  atmp(1:dim, 1:dim) = 0.d+0
  do i = 1, dim
     do j = 1, dim
        do k = 1, dim
           do l = 1, dim
              atmp(i, j) = atmp(i, j) + uvec(k, i) * amat(k, l) * vvec(l, j)
           end do
        end do
     end do
  end do

!debug  write(6, "('util_youla_real: a = Ut * A * V <=> A = U * a * Vt')")
!debug  do i = 1, dim
!debug     do j = 1, dim
!debug        write(6, "(f20.10)", advance = 'no') atmp(i, j)
!debug     end do
!debug     write(6, *)
!debug  end do

  do ii = 1, dim2
     i = 2 * ii - 1
     aeig(ii) = atmp(i, i+1)
  end do

  write(6, "('util_youla_real. eigenvalues:')")
  do ii = 1, dim2
     write(6, "(i10, f20.10)") ii, aeig(ii)
  end do

  deallocate(atmp)

end subroutine util_youla_real
!################################################################################
