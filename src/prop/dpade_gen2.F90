!######################################################################
subroutine dpade_gen2(coeff0, coeff1, tinv, tpiv)

! Generate LU factorization of coeff0 + coeff1 * T
! T is kinetic + nucleus-electron operator.

  use, intrinsic :: iso_c_binding
  use mod_const, only : half, czero, runit, iunit
  use mod_control, only : fedvr_normalized
  use mod_rad, only : nrad, wrad, ndvr
  use mod_sph, only : lmax1
  use mod_bas, only : tmat

  implicit none
  complex(c_double_complex), intent(in) :: coeff0, coeff1
  complex(c_double_complex), intent(out) :: tinv(1:(3*ndvr+1), 1:(nrad-1), 0:lmax1)
  integer(c_int), intent(out) :: tpiv(1:(nrad-1), 0:lmax1)

  integer(c_int) :: irad, jrad, l, jll, jul, jb1, jb2, dim, ld2, info
  complex(c_double_complex) :: fac

  dim = nrad - 1
  ld2 = 3 * ndvr + 1
  tinv(1:ld2, 1:dim, 0:lmax1) = czero

  do l = 0, lmax1
     do irad = 1, nrad - 1
        tinv(2*ndvr+1, irad, l) = coeff0
        jll = max(1,        irad - ndvr)
        jul = min(nrad - 1, irad + ndvr)
        if (fedvr_normalized) then
           do jrad = jll, jul
              jb1 =     ndvr + 1 + jrad - irad
              jb2 = 2 * ndvr + 1 + jrad - irad
              tinv(jb2, irad, l) = tinv(jb2, irad, l) + tmat(jb1, irad, l) * coeff1
           end do
        else
           do jrad = jll, jul
              jb1 =     ndvr + 1 + jrad - irad
              jb2 = 2 * ndvr + 1 + jrad - irad
              tinv(jb2, irad, l) = tinv(jb2, irad, l) + tmat(jb1, irad, l) * coeff1 / wrad(jrad)
           end do
        end if
     end do

     call zgbtrf(dim, dim, ndvr, ndvr, tinv(1,1,l), ld2, tpiv(1,l), info)
     if (info .ne. 0) then
        write(6, "('dpade_gen2-zgbtrf: l = ', i5, ' info = ', i5)") l, info
        stop 'Error of zgbtrf.'
     end if
  end do

end subroutine dpade_gen2
!######################################################################
