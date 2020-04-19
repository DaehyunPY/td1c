!######################################################################
subroutine bas_gen_d2ll(d2ll)

  use, intrinsic :: iso_c_binding
  use mod_sph, only : lmax2
  use mod_const, only : zero, one, two
  use mod_control, only : fedvr_normalized
  use mod_rad, only : nrad, ndvr, xrad, wrad, radk

  implicit none
  real(c_double), intent(out) :: d2ll(1:(ndvr+1), 1:(nrad-1), 0:lmax2)

  integer(c_int) :: ifun, irad, jrad, krad, l, l2, m
  integer(c_int) :: jll, jul, jul2, jb1, jb2
  integer(c_int) :: dim, ld, info
  real(c_double) :: radk2, angk, angk2
!debug
!  integer(c_int) :: ldz = 1
!  real(c_double), allocatable :: d2eig(:), work(:), zz(:,:)
!debug

  dim = nrad - 1
  ld = ndvr + 1

  d2ll(1:ld, 1:dim, 0:lmax2) = zero
  if (fedvr_normalized) then
     do l = 0, lmax2
        l2 = l * (l + 1)
        do irad = 1, nrad - 1
           angk = one / xrad(irad)
           angk2 = angk * angk * dble(l2);
           d2ll(1, irad, l) = angk2
        end do
     end do
  else
     do l = 0, lmax2
        l2 = l * (l + 1)
        do irad = 1, nrad - 1
           angk = one / xrad(irad)
           angk2 = angk * angk * dble(l2);
           d2ll(1, irad, l) = angk2 * wrad(irad)
        end do
     end do
  end if

  do l = 0, lmax2
     l2 = l * (l + 1)
     do irad = 1, nrad - 1
        jll = max(1,   irad)
        jul = min(dim, irad + ndvr)
        do jrad = jll, jul
           jb1 = 1 + jrad - irad + ndvr
           jb2 = 1 + jrad - irad
           radk2 = radk(jb1, irad) * two
           d2ll(jb2, irad, l) = d2ll(jb2, irad, l) + radk2
        end do
     end do
  end do

!debug
!  do l = 0, lmax2
!     write(6, "('bas_gen_d2ll: Laplacian, l = ', i10)") l
!     do irad = 1, nrad - 1
!        jll  = max(1,   irad - ndvr)
!        jul  = min(dim, irad)
!        jul2 = min(dim, irad + ndvr)
!        do jrad = 1, nrad - 1
!           if (jrad >= jll .and. jrad <= jul) then
!              ib2 = 1 + irad - jrad
!              write(6, "(f15.10)", advance = 'no') dble(d2ll(ib2, jrad, l))
!           else if (jrad >= jll .and. jrad <= jul2) then
!              jb2 = 1 + jrad - irad
!              write(6, "(f15.10)", advance = 'no') dble(d2ll(jb2, irad, l))
!           else
!              write(6, "('      xxx      ')", advance = 'no')
!           end if
!        end do
!        write(6, *)
!     end do
!  end do     
!debug

!  ldz = dim
!  allocate(d2eig(1:dim))
!  allocate(work(1:(3*dim-2)))
!  allocate(zz(1:ldz, 1:dim))
!  do l = 0, lmax2
!     call dsbev('V', 'L', dim, ndvr, d2ll(1,1,l), ld, d2eig, zz, ldz, work, info)
!     write(6, "('eigenvalus of d2ll: ', 2i5)") l, info
!     do irad = 1, dim
!        write(6, "(i5, f20.10)") irad, d2eig(irad)
!     end do
!!     write(6, "('first eigenvectors of d2ll: ')")
!!     do irad = 1, dim
!!        write(6, "(i5, 10f10.5)") irad, zz(irad, 1:10)
!!     end do
!  end do
!  deallocate(zz)
!  deallocate(work)
!  deallocate(d2eig)
!  stop
!debug

  do l = 0, lmax2
!DEBUG
!     write(6, "('d2ll: l    = ', i5)") l
!     write(6, "('d2ll: dim  = ', i5)") dim
!     write(6, "('d2ll: ndvr = ', i5)") ndvr
!     write(6, "('d2ll: ld   = ', i5)") ld
!     do irad = 1, nrad-1
!        write(6,"('d2ll:',i5)",advance='no') irad
!        do m = 1, ndvr+1
!           write(6,"(f10.5)",advance='no') d2ll(m,irad,l)
!        end do
!        write(6,*)
!     end do
!     stop
!DEBUG
     call dpbtrf('L', dim, ndvr, d2ll(1,1,l), ld, info)
     if (info /= 0) then
        write(6, "('dpbtrf: l = ', i5, ' info = ', i15)") l, info
        stop
     end if
!     write(6, "('bas_gen_d2ll: LL-factorized Laplacian, l = ', i10)") l
!     do irad = 1, nrad - 1
!        jll = max(1, irad - ndvr)
!        jul = min(dim, irad)
!        do jrad = 1, nrad - 1
!           if (jrad >= jll .and. jrad <= jul) then
!              ib = 1 + irad - jrad
!              write(6, "(f15.10)", advance = 'no') dble(d2ll(ib, jrad, l))
!           else
!              write(6, "('      xxx      ')", advance = 'no')
!           end if
!        end do
!        write(6, *)
!     end do
  end do

end subroutine bas_gen_d2ll
!######################################################################
