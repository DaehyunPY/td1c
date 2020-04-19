!2016/10/22 Yuki Orimo Changed
!######################################################################
subroutine bas_gen_kmat(kmat)
!
  use, intrinsic :: iso_c_binding
  use mod_sph, only : lmax1, mmax1
  use mod_const, only : zero, one, half
  use mod_control, only : fedvr_normalized
  use mod_rad, only : nrad, ndvr, xrad, wrad, radk
!
! atomic hamiltonian in banded storage
!
! example: KL = 2, KU = 1:
! <--------- input ---------->     <--------- output --------->
!  *    *    *    +    +    +       *    *    *   u14  u25  u36
!  *    *    +    +    +    +       *    *   u13  u24  u35  u46
!  *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
! a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
! a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *
! a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *
!
  implicit none
  complex(c_double_complex), intent(out) :: kmat(1:(2*ndvr+1), 1:(nrad-1), 0:lmax1)
  integer(c_long) :: ifun, irad, jrad, krad, l, l2, m, jll, jul, jb1
  integer(c_long) :: dim, nsub, ld1
  real(c_double) :: oor, zor, zor1

  dim = nrad - 1
  nsub = ndvr
  ld1 = 2 * nsub + 1

  kmat (1:ld1, 1:dim, 0:lmax1) = zero
  if (fedvr_normalized) then
     do l = 0, lmax1
        l2 = l * (l + 1)
        do irad = 1, nrad - 1
           oor = one / xrad(irad)
           zor = oor * dble(l2) * oor * half;
           kmat (1 + nsub, irad, l) = zor
        end do
     end do
  else
     do l = 0, lmax1
        l2 = l * (l + 1)
        do irad = 1, nrad - 1
           oor = one / xrad(irad)
           zor = oor * dble(l2) * oor * half;
           kmat (1 + nsub, irad, l) = zor * wrad(irad)
        end do
     end do
  end if

  do l = 0, lmax1
     l2 = l * (l + 1)
     do irad = 1, nrad - 1
        jll = max(1,        irad - nsub)
        jul = min(nrad - 1, irad + nsub)
        do jrad = jll, jul
           jb1 =  nsub + 1 + jrad - irad
           kmat(jb1, irad, l) = kmat(jb1, irad, l) + radk(jb1, irad)
        end do
     end do
  end do

  !DEBUG
  !DEBUG call bas_test_keig(nrad, lmax1, mmax1, ndvr, xrad, kmat)
  !DEBUG

end subroutine bas_gen_kmat
!######################################################################
subroutine bas_test_keig(nrad, lmax1, mmax1, ndvr, xrad, kmat)
!
  use, intrinsic :: iso_c_binding

! atomic hamiltonian and its lu decomposition in banded matrices
!
  implicit none
  integer(c_long), intent(in) :: nrad, lmax1, mmax1, ndvr
  real(c_double), intent(in) :: xrad(0:nrad)
  real(c_double), intent(inout) :: kmat(1:(ndvr+1), 1:(nrad-1), 0:lmax1)

  integer(c_long) :: ifun, irad, jrad, krad, l, l2, m, jll, jul, jul2, ib, jb
  integer(c_long) :: dim, nsub, ld, info
  real(c_double) :: kang, trad2, kang2
  integer(c_long) :: ldz = 1
  real(c_double), allocatable :: keig(:), work(:), zz(:,:)

  dim = nrad - 1
  nsub = ndvr
  ld = nsub + 1
  ldz = dim

  allocate(keig(1:dim))
  allocate(work(1:(3*dim-2)))
  allocate(zz(1:ldz, 1:dim))
  do l = 0, 0
     call dsbev('V', 'L', dim, nsub, kmat(1,1,l), ld, keig, zz, ldz, work, info)
     write(6, "('eigenvalus of kmat: ', 2i5)") l, info
     do irad = 1, dim
        write(6, "(i5, f20.10)") irad, keig(irad)
     end do
!     write(6, "('first eigenvectors of kout: ')")
!     do irad = 1, dim
!        write(6, "(i5, 10f10.5)") irad, zz(irad, 1:10)
!     end do
  end do
  deallocate(zz)
  deallocate(work)
  deallocate(keig)
  stop "STOP for debug @ bas_gen_keig."

end subroutine bas_test_keig
!######################################################################
