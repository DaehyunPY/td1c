!######################################################################
subroutine adi_gen_tadi(icomp, dokin, donuc, dtime, tadi1, tadi2, tpiv2)

  use, intrinsic :: iso_c_binding
  use mod_const, only : one, half, czero, runit, iunit
  use mod_rad, only : nrad, xrad, radk, ndvr
  use mod_sph, only : lmax1
  use mod_bas, only : znuc

!
! atomic hamiltonian and its lu decomposition in banded matrices
!
  implicit none
  integer(c_long), intent(in) :: icomp, dokin, donuc
  real(c_double), intent(in) :: dtime
  complex(c_double_complex), intent(out) :: tadi1(1:(2*ndvr+1), 1:(nrad-1), 0:lmax1)
  complex(c_double_complex), intent(out) :: tadi2(1:(3*ndvr+1), 1:(nrad-1), 0:lmax1)
  integer(c_long), intent(out) :: tpiv2(1:(nrad-1), 0:lmax1)

  integer(c_long) :: irad, jrad, l, l2, jll, jul, ib1, ib2
  integer(c_long) :: dim, nsub, ld1, ld2, info
  complex(c_double_complex) :: fac, fradk
  real(c_double) :: oor, zor

!debug
!  integer(c_long) :: ld
!  real(c_double), allocatable :: teig(:), rwork(:)  
!  complex(c_double_complex), allocatable :: f1(:,:,:), work(:), zz(:,:)
!debug

  dim = nrad - 1
  nsub = ndvr
  ld1 = 2 * nsub + 1
  ld2 = 3 * nsub + 1
  if (icomp == 1) then
     fac = iunit * dtime * half
  else
     fac = runit * dtime * half
  end if

  tadi1(1:ld1, 1:dim, 0:lmax1) = czero
  tadi2(1:ld2, 1:dim, 0:lmax1) = czero
  do l = 0, lmax1
     l2 = l * (l + 1)
     do irad = 1, nrad - 1
        ib1 = 1 + nsub
        ib2 = 1 + nsub * 2
        if (donuc .ne. 0) then
           oor = one / xrad(irad)
           zor = oor * (-dble(znuc) + dble(l2) * oor * half);
        else
           oor = one / xrad(irad)
           zor = oor * dble(l2) * oor * half;
        end if
        tadi1(ib1, irad, l) = runit - zor * fac
        tadi2(ib2, irad, l) = runit + zor * fac

!bug        if (donuc .ne. 0) then
!bug           oor = 1.0d+0 / xrad(irad)
!bug           zor = oor * (-dble(znuc) + dble(l2) * oor * 0.5d+0);
!bug           !reg zor1 = oor * (0.d+0 + dble(l2) * oor * 0.5d+0);
!bug           tadi1(ib1, irad, l) = runit - zor * fac
!bug           tadi2(ib2, irad, l) = runit + zor * fac
!bug        else
!bug           tadi1(ib1, irad, l) = runit
!bug           tadi2(ib2, irad, l) = runit
!bug        end if

        if (dokin .ne. 0) then
           jll = max(1,        irad - nsub)
           jul = min(nrad - 1, irad + nsub)
           do jrad = jll, jul
              ib1 =     nsub + 1 + irad - jrad
              ib2 = 2 * nsub + 1 + irad - jrad

              fradk = fac * radk(ib1, jrad)
              tadi1(ib1, jrad, l) = tadi1(ib1, jrad, l) - fradk
              tadi2(ib2, jrad, l) = tadi2(ib2, jrad, l) + fradk
           end do
        end if
     end do
  end do

!?? same as above ??  tadi1(1:ld1, 1:dim, 0:lmax1) = czero
!?? same as above ??  tadi2(1:ld2, 1:dim, 0:lmax1) = czero
!?? same as above ??  do l = 0, lmax1
!?? same as above ??     l2 = l * (l + 1)
!?? same as above ??     do irad = 1, nrad - 1
!?? same as above ??        oor = 1.0d+0 / xrad(irad)
!?? same as above ??        zor = oor * (-dble(znuc) + dble(l2) * oor * 0.5d+0);
!?? same as above ??
!?? same as above ??        ib1 = 1 + nsub
!?? same as above ??        ib2 = 1 + nsub * 2
!?? same as above ??        tadi1(ib1, irad, l) = runit - zor * fac
!?? same as above ??        tadi2(ib2, irad, l) = runit + zor * fac
!?? same as above ??
!?? same as above ??        jll = max(1, irad - nsub)
!?? same as above ??        jul = min(nrad - 1, irad + nsub)
!?? same as above ??        do jrad = jll, jul
!?? same as above ??           ib1 =     nsub + 1 + irad - jrad
!?? same as above ??           ib2 = 2 * nsub + 1 + irad - jrad
!?? same as above ??           ij = (2 * nsub + 1) * irad + (jrad - (irad - nsub)) + 1 - nsub
!?? same as above ??           tadi1(ib1, jrad, l) = tadi1(ib1, jrad, l) - kmat(ij) * fac
!?? same as above ??           tadi2(ib2, jrad, l) = tadi2(ib2, jrad, l) + kmat(ij) * fac
!?? same as above ??        end do
!?? same as above ??     end do
!?? same as above ??  end do

!debug  do l = 0, lmax1
!debug     write(6, "('adi_gen_tadi: 1 - i*T*dtime/2, l = ', i10, ' (real)')") l
!debug     do irad = 1, nrad - 1
!debug        jll = max(1, irad - nsub)
!debug        jul = min(dim, irad + nsub)
!debug        do jrad = 1, nrad - 1
!debug           if (jrad >= jll .and. jrad <= jul) then
!debug              ib1 = nsub + 1 + irad - jrad
!debug              write(6, "(f15.10)", advance = 'no') dble(tadi1(ib1, jrad, l))
!debug           else
!debug              write(6, "('      xxx      ')", advance = 'no')
!debug           end if
!debug        end do
!debug        write(6, *)
!debug     end do
!debug     write(6, "('adi_gen_tadi: 1 - i*T*dtime/2, l = ', i10, ' (imag)')") l
!debug     do irad = 1, nrad - 1
!debug        jll = max(1, irad - nsub)
!debug        jul = min(dim, irad + nsub)
!debug        do jrad = 1, nrad - 1
!debug           if (jrad >= jll .and. jrad <= jul) then
!debug              ib1 = nsub + 1 + irad - jrad
!debug              write(6, "(f15.10)", advance = 'no') aimag(tadi1(ib1, jrad, l))
!debug           else
!debug              write(6, "('      xxx      ')", advance = 'no')
!debug           end if
!debug        end do
!debug        write(6, *)
!debug     end do
!debug  end do
     

  do l = 0, lmax1
     call zgbtrf(dim, dim, nsub, nsub, tadi2(1,1,l), ld2, tpiv2(1,l), info)
     if (info .ne. 0) then
        write(6, "('adi_gen_tadi-zgbtrf: l = ', i5, ' info = ', i5)") l, info
        stop 'Error of zgbtrf.'
     end if
!test     call zpotrf('L', dim, tadi2(1,1,l), ld2, info)
!test     write(6, "('zpotrf: info = ', i5)") info
!     write(6, "('adi_gen_tadi: lu facts of 1 + i*T*dtime/2, l = ', i10, ' (real)')") l
!     do irad = 1, nrad - 1
!        jll = max(1, irad - nsub)
!        jul = min(dim, irad + nsub)
!        do jrad = 1, nrad - 1
!           if (jrad >= jll .and. jrad <= jul) then
!              ib2 = 2 * nsub + 1 + irad - jrad
!              write(6, "(f15.10)", advance = 'no') dble(tadi2(ib2, jrad, l))
!           else
!              write(6, "('      xxx      ')", advance = 'no')
!           end if
!        end do
!        write(6, *)
!     end do
!     write(6, "('adi_gen_tadi: lu facts of 1 + i*T*dtime/2, l = ', i10, ' (imag)')") l
!     do irad = 1, nrad - 1
!        jll = max(1, irad - nsub)
!        jul = min(dim, irad + nsub)
!        do jrad = 1, nrad - 1
!           if (jrad >= jll .and. jrad <= jul) then
!              ib2 = 2 * nsub + 1 + irad - jrad
!              write(6, "(f15.10)", advance = 'no') aimag(tadi2(ib2, jrad, l))
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
!  allocate(f1(1:(ndvr+1), 1:dim, 0:lmax1))
!  f1(1:(ndvr+1), 1:dim, 0:lmax1) = czero
!  do l = 0, lmax1
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
!  do l = 0, lmax1
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


end subroutine adi_gen_tadi
!######################################################################
!   example: KL = 2, KU = 1:
!   <--------- input ---------->     <--------- output --------->
!    *    *    *    +    +    +       *    *    *   u14  u25  u36
!    *    *    +    +    +    +       *    *   u13  u24  u35  u46
!    *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
!   a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
!   a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *
!   a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *
