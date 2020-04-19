!######################################################################
subroutine pade_gen2(icomp, dokin, donuc, nrad, nsph, lmax, mmax, noffd, &
     znuc, dtime, alpha, xrad, kmat, tinv, tpiv)
  use, intrinsic :: iso_c_binding

!
! Generate LU factorization of 
! -i*T*dt/2 - alpha (icomp == 1), or
! -1*T*dt/2 - alpha (icomp == 0), where
! T is kinetic + nucleus-electron operator.
!
  implicit none
  integer(c_int), intent(in) :: icomp, dokin, donuc
  integer(c_int), intent(in) :: nrad, nsph, lmax, mmax, noffd, znuc
  real(c_double), intent(in) :: dtime
  complex(c_double_complex), intent(in) :: alpha
  real(c_double), intent(in) :: xrad(0:nrad)
  real(c_double), intent(in) :: kmat(1:*)
  complex(c_double_complex), intent(out) :: tinv(1:(3*noffd+1), 1:(nrad-1), 0:lmax)
  integer(c_int), intent(out) :: tpiv(1:(nrad-1), 0:lmax)

  integer(c_int) :: ifun, irad, jrad, krad, l, l2, m, lm, jll, jul, ij, mlim
  integer(c_int) :: ib, ib1, ib2

  integer(c_int) :: dim, nsub, ld1, ld2, info
  real(c_double) :: oor, zor, zor1
  complex(c_double_complex) :: fac
  complex(c_double_complex), parameter :: czero = (0.d+0, 0.d+0)
  complex(c_double_complex), parameter :: runit = (1.d+0, 0.d+0)
  complex(c_double_complex), parameter :: iunit = (0.d+0, 1.d+0)

!debug
!  integer(c_int) :: ld
!  real(c_double), allocatable :: teig(:), rwork(:)  
!  complex(c_double_complex), allocatable :: f1(:,:,:), work(:), zz(:,:)
!debug


  dim = nrad - 1
  nsub = noffd
  ld1 = 2 * nsub + 1
  ld2 = 3 * nsub + 1
  if (icomp == 1) then
     fac = iunit * dtime / alpha
  else
     fac = runit * dtime / alpha
  end if

  tinv(1:ld2, 1:dim, 0:lmax) = czero
  do l = 0, lmax
     l2 = l * (l + 1)
     do irad = 1, nrad - 1
        ib1 = 1 + nsub
        ib2 = 1 + nsub * 2
        if (donuc .ne. 0) then
           oor = 1.0d+0 / xrad(irad)
           zor = oor * (-dble(znuc) + dble(l2) * oor * 0.5d+0);
        else
           oor = 1.0d+0 / xrad(irad)
           zor = oor * dble(l2) * oor * 0.5d+0;
        end if
        tinv(ib2, irad, l) = runit + zor * fac

        if (dokin .ne. 0) then
           jll = max(1, irad - nsub)
           jul = min(nrad - 1, irad + nsub)
           do jrad = jll, jul
              ib1 =     nsub + 1 + irad - jrad
              ib2 = 2 * nsub + 1 + irad - jrad
              ij = (2 * nsub + 1) * irad + (jrad - (irad - nsub)) + 1 - nsub
              tinv(ib2, jrad, l) = tinv(ib2, jrad, l) + kmat(ij) * fac
           end do
        end if
     end do
  end do

  do l = 0, lmax
     call zgbtrf(dim, dim, nsub, nsub, tinv(1,1,l), ld2, tpiv(1,l), info)
     if (info .ne. 0) then
        write(6, "('pade_gen2-zgbtrf: l = ', i5, ' info = ', i5)") l, info
        stop 'Error of zgbtrf.'
     end if
!test     call zpotrf('L', dim, tinv(1,1,l), ld2, info)
!test     write(6, "('zpotrf: info = ', i5)") info
!     write(6, "('pade_gen2: lu facts of 1 + i*T*dtime/2, l = ', i10, ' (real)')") l
!     do irad = 1, nrad - 1
!        jll = max(1, irad - nsub)
!        jul = min(dim, irad + nsub)
!        do jrad = 1, nrad - 1
!           if (jrad >= jll .and. jrad <= jul) then
!              ib2 = 2 * nsub + 1 + irad - jrad
!              write(6, "(f15.10)", advance = 'no') dble(tinv(ib2, jrad, l))
!           else
!              write(6, "('      xxx      ')", advance = 'no')
!           end if
!        end do
!        write(6, *)
!     end do
!     write(6, "('pade_gen2: lu facts of 1 + i*T*dtime/2, l = ', i10, ' (imag)')") l
!     do irad = 1, nrad - 1
!        jll = max(1, irad - nsub)
!        jul = min(dim, irad + nsub)
!        do jrad = 1, nrad - 1
!           if (jrad >= jll .and. jrad <= jul) then
!              ib2 = 2 * nsub + 1 + irad - jrad
!              write(6, "(f15.10)", advance = 'no') aimag(tinv(ib2, jrad, l))
!           else
!              write(6, "('      xxx      ')", advance = 'no')
!           end if
!        end do
!        write(6, *)
!     end do
  end do

!debug
!  allocate(teig(1:dim))
!  allocate(work(1:dim))
!  allocate(rwork(1:(max(1,3*dim-2))))
!  allocate(zz(1:dim, 1:dim))
!  allocate(f1(1:(noffd+1), 1:dim, 0:lmax))
!  f1(1:(noffd+1), 1:dim, 0:lmax) = czero
!  do l = 0, lmax
!     l2 = l * (l + 1)
!     do irad = 1, nrad - 1
!        oor = 1.0d+0 / xrad(irad)
!        zor = oor * (-dble(znuc) + dble(l2) * oor * 0.5d+0);
!
!        ib1 = 1
!        f1(ib1, irad, l) = zor
!
!        jll = max(1, irad - nsub)
!        jul = min(dim, irad)
!        do jrad = jll, jul
!           ib1 = 1 + irad - jrad
!           ij = (2 * nsub + 1) * irad + (jrad - (irad - nsub)) + 1 - nsub
!           f1(ib1, jrad, l) = f1(ib1, jrad, l) + kmat(ij)
!        end do
!     end do
!  end do
!
!  ld = nsub + 1
!  do l = 0, lmax
!     call zhbev('V', 'L', dim, nsub, f1(1,1,l), ld, &
!          & teig, zz, dim, work, rwork, info)
!
!     write(6, "('eigenvalues of tmat: ', 2i5)") l, info
!     do irad = 1, dim
!        write(6, "(i5, f20.10)") irad, teig(irad)
!     end do
!!     write(6, "('first eigenvectors of tmat: ')")
!!     do irad = 1, dim
!!        write(6, "(i5, 10f10.5)") irad, zz(irad, 1:10)
!!     end do
!  end do
!  deallocate(f1)
!  deallocate(zz)
!  deallocate(rwork)
!  deallocate(work)
!  deallocate(teig)
!!  stop
!debug


end subroutine pade_gen2
!######################################################################
!   example: KL = 2, KU = 1:
!   <--------- input ---------->     <--------- output --------->
!    *    *    *    +    +    +       *    *    *   u14  u25  u36
!    *    *    +    +    +    +       *    *   u13  u24  u35  u46
!    *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
!   a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
!   a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *
!   a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *
