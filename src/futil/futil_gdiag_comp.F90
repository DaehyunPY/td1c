!################################################################################
subroutine futil_gdiag_comp(dsc, n, a, uleft, uright)

  use, intrinsic :: iso_c_binding

  implicit none
  logical, intent(in) :: dsc
  integer(c_long), intent(in) :: n
  complex(c_double_complex), intent(inout) :: a(1:n, 1:n)
  complex(c_double_complex), intent(out) :: uleft(1:n, 1:n)
  complex(c_double_complex), intent(out) :: uright(1:n, 1:n)

  integer(c_long) :: info, lwork, i, j, k, len
  complex(c_double_complex) :: clwork, test

  integer(c_long), allocatable :: index(:)
  complex(c_double_complex), allocatable :: tmpa(:, :)
  complex(c_double_complex), allocatable :: eig(:), work(:), vecl(:,:), vecr(:,:), plr(:)
  real(c_double), allocatable :: rwork(:), deig(:)

  len = 2 * n
  allocate(index(n))
  allocate(rwork(len))
  allocate(tmpa(n, n))
  allocate(vecl(n, n))
  allocate(vecr(n, n))
  allocate(eig(n))
  allocate(deig(n))
  allocate(plr(n))
  tmpa(1:n, 1:n) = a(1:n, 1:n)

!debug
!debug  write(6, "('futil_gdiag_comp: n = ', i10)") n
!debug  write(6, "('A: A')")
!debug  write(6, "(2f20.10)") tmpa(1:n, 1:n)
!debug
  call zgeev("v", "v", n, tmpa, n, eig, vecl, n, vecr, n, clwork, -1, rwork, info)
  lwork = int(clwork)

  allocate(work(lwork))
  call zgeev("v", "v", n, tmpa, n, eig, vecl, n, vecr, n, work, lwork, rwork, info)

  if (info /= 0) then
     write(6, "('futil_gdiag_comp: info = ', i20)") info
     stop 'error in futil_gdiag_comp.'
  else
     uleft (1:n, 1:n) = (0.d+0, 0.d+0)
     uright(1:n, 1:n) = (0.d+0, 0.d+0)
     a     (1:n, 1:n) = (0.d+0, 0.d+0)
!     if (dsc) then
!        do i = 1, n
!           a(i, i) = eig(n - i + 1)
!           uleft (1:n, i) = vecl(1:n, n - i + 1)
!           uright(1:n, i) = vecr(1:n, n - i + 1)
!        end do
!     else
!        do i = 1, n
!           a(i, i) = eig(i)
!           uleft (1:n, i) = vecl(1:n, i)
!           uright(1:n, i) = vecr(1:n, i)
!        end do
!     end if

     deig(1:n) = dble(eig(1:n))
     call futil_index_double(dsc, n, deig, index)
     do i = 1, n
        j = index(i)
        a(i, i) = eig(j)
        uleft (1:n, i) = vecl(1:n, j)
        uright(1:n, i) = vecr(1:n, j)
     end do
  end if

  ! renormalize eigenvector to hold biorthonormality
  ! write(6, "('gdiag_comp: biorthogonality check')")
  do i = 1, n
     do j = 1, n
        test = (0.d+0, 0.d+0)
        do k = 1, n
           test = test + conjg(uleft(k, i)) * uright(k, j)
        end do
        if (i == j) plr(i) = test
        ! write(6, "(2i5, 2f20.10)") i, j, test
     end do
  end do
  do i = 1, n
     do j = 1, n
        test = sqrt(plr(j))
        uright(i, j) = uright(i, j) / test
        uleft(i, j) = uleft(i, j) / conjg(test)
     end do
  end do

  deallocate(work)
  deallocate(plr)
  deallocate(deig)
  deallocate(eig)
  deallocate(vecr)
  deallocate(vecl)
  deallocate(tmpa)
  deallocate(rwork)
  deallocate(index)

end subroutine futil_gdiag_comp
!################################################################################
