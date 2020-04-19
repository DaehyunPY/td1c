!######################################################################
subroutine h1rat_gen(dimd, cd0, cd1, tinv, tpiv)

! Generate LU factorization of 
! cd0[i] + cd1[i] * T, where
! T is kinetic + nucleus-electron operator, for
! i = 1, 2, , ..., dimd.

  use, intrinsic :: iso_c_binding
  use mod_bas, only : tmat
  use mod_sph, only : lmax1
  use mod_const, only : czero
  use mod_rad, only : nrad, wrad, ndvr
  use mod_control, only : fedvr_normalized, h1rat_thresh

  implicit none
  integer(c_int), intent(in) :: dimd
  complex(c_double_complex), intent(in) :: cd0(1:dimd)
  complex(c_double_complex), intent(in) :: cd1(1:dimd)
  complex(c_double_complex), intent(out) :: tinv(1:(3*ndvr+1), 1:(nrad-1), 0:lmax1, 1:dimd)
  integer(c_int), intent(out) :: tpiv(1:(nrad-1), 0:lmax1, 1:dimd)

  integer(c_int) :: idim, irad, jrad, l, jll, jul, jb1, jb2, dim, ld2, info
  complex(c_double_complex) :: fac

  dim = nrad - 1
  ld2 = 3 * ndvr + 1
  tinv(1:ld2, 1:dim, 0:lmax1, 1:dimd) = czero

  do idim = 1, dimd
     do l = 0, lmax1
        do irad = 1, nrad - 1
           tinv(2*ndvr+1, irad, l, idim) = cd0(idim)
           jll = max(1,        irad - ndvr)
           jul = min(nrad - 1, irad + ndvr)
           if (fedvr_normalized) then
              do jrad = jll, jul
                 jb1 =     ndvr + 1 + jrad - irad
                 jb2 = 2 * ndvr + 1 + jrad - irad
                 tinv(jb2, irad, l, idim) = tinv(jb2, irad, l, idim) + tmat(jb1, irad, l) * cd1(idim)
              end do
           else
              stop 'h1rat_gen: unnormalized basis not supported.'
           end if
        end do

        if (abs(cd1(idim)) > h1rat_thresh) then
           call zgbtrf(dim, dim, ndvr, ndvr, tinv(1,1,l,idim), ld2, tpiv(1,l,idim), info)
           if (info .ne. 0) then
              write(6, "('h1rat_gen-zgbtrf: idim = ', i5, ' l = ', i5, ' info = ', i5)") idim, l, info
              stop 'Error of zgbtrf.'
           end if
        end if
     end do
  end do

end subroutine h1rat_gen
!######################################################################
