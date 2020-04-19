!######################################################################
subroutine h1itr_gen_cninv(icomp, dtime, cnpiv, cninv)

! Generate LU factorization of 
! 1+i*T*dt/2 (icomp == 1), or
! 1+1*T*dt/2 (icomp == 0), where
! T is kinetic + nucleus-electron operator.

  use, intrinsic :: iso_c_binding
  use mod_const, only : half, czero, runit, iunit
  use mod_control, only : fedvr_normalized
  use mod_rad, only : nrad, wrad, ndvr
  use mod_sph, only : lmax1
  use mod_bas, only : tmat

  implicit none
  integer(c_int), intent(in) :: icomp
  real(c_double), intent(in) :: dtime
  integer(c_int), intent(out) :: cnpiv(1:(nrad-1), 0:lmax1)
  complex(c_double_complex), intent(out) :: cninv(1:(3*ndvr+1), 1:(nrad-1), 0:lmax1)

  integer(c_int) :: irad, jrad, l, jll, jul, jb1, jb2, dim, ld2, info
  complex(c_double_complex) :: fac

  !DEBUG
  !call h1itr_gen_cninv2(icomp, dtime, cnpiv, cninv)
  !return
  !DEBUG

  dim = nrad - 1
  ld2 = 3 * ndvr + 1
  cninv(1:ld2, 1:dim, 0:lmax1) = czero

  if (icomp == 1) then
     fac = iunit * dtime * half
  else
     fac = runit * dtime * half
  end if

  do l = 0, lmax1
     do irad = 1, nrad - 1
        if (fedvr_normalized) then
           cninv(2*ndvr+1, irad, l) = runit
        else
           cninv(2*ndvr+1, irad, l) = wrad(irad)
        end if
        jll = max(1,        irad - ndvr)
        jul = min(nrad - 1, irad + ndvr)
        do jrad = jll, jul
           jb1 =     ndvr + 1 + jrad - irad
           jb2 = 2 * ndvr + 1 + jrad - irad
           cninv(jb2, irad, l) = cninv(jb2, irad, l) + fac * tmat(jb1, irad, l)
        end do
     end do

     call zgbtrf(dim, dim, ndvr, ndvr, cninv(1,1,l), ld2, cnpiv(1,l), info)
     if (info .ne. 0) then
        write(6, "('h1itr_gen_cninv-zgbtrf: l = ', i5, ' info = ', i5)") l, info
        stop 'Error of zgbtrf.'
     end if
  end do

end subroutine h1itr_gen_cninv
!######################################################################
subroutine h1itr_gen_cninv2(icomp, dtime, cnpiv, cninv)

! Generate LU factorization of 
! 1+i*T*dt/2 (icomp == 1), or
! 1+1*T*dt/2 (icomp == 0), where
! T is kinetic + nucleus-electron operator.

  use, intrinsic :: iso_c_binding
  use mod_const, only : half, czero, runit, iunit
  use mod_control, only : fedvr_normalized
  use mod_rad, only : nrad, wrad, ndvr
  use mod_sph, only : lmax1
  use mod_bas, only : tmat

  implicit none
  integer(c_int), intent(in) :: icomp
  real(c_double), intent(in) :: dtime
  integer(c_int), intent(out) :: cnpiv(1:(nrad-1), 0:lmax1)
  complex(c_double_complex), intent(out) :: cninv(1:(3*ndvr+1), 1:(nrad-1), 0:lmax1)


  integer(c_int) :: irad, jrad, l, jll, jul, jb1, jb2, dim, ld2, info
  complex(c_double_complex) :: fac

  dim = nrad - 1
  ld2 = 3 * ndvr + 1
  cninv(1:ld2, 1:dim, 0:lmax1) = czero

  if (icomp == 1) then
     fac = iunit * dtime * half
  else
     fac = runit * dtime * half
  end if

  do l = 0, lmax1
     do irad = 1, nrad - 1
        cninv(2*ndvr+1, irad, l) = runit
        jll = max(1,        irad - ndvr)
        jul = min(nrad - 1, irad + ndvr)
        if (fedvr_normalized) then
           do jrad = jll, jul
              jb1 =     ndvr + 1 + jrad - irad
              jb2 = 2 * ndvr + 1 + jrad - irad
              cninv(jb2, irad, l) = cninv(jb2, irad, l) + fac * tmat(jb1, irad, l)
           end do
        else
           do jrad = jll, jul
              jb1 =     ndvr + 1 + jrad - irad
              jb2 = 2 * ndvr + 1 + jrad - irad
              cninv(jb2, irad, l) = cninv(jb2, irad, l) + fac * tmat(jb1, irad, l) / wrad(irad)
           end do
        end if
     end do

     call zgbtrf(dim, dim, ndvr, ndvr, cninv(1,1,l), ld2, cnpiv(1,l), info)
     if (info .ne. 0) then
        write(6, "('h1itr_gen_cninv-zgbtrf: l = ', i5, ' info = ', i5)") l, info
        stop 'Error of zgbtrf.'
     end if
  end do

end subroutine h1itr_gen_cninv2
!######################################################################
