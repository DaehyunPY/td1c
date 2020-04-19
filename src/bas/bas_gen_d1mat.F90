!######################################################################
subroutine bas_gen_d1mat(d1mat)
!
  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad, ndvr, xrad, radp, ecs_flag, irad_ecs, theta, wrad, cwrad, mapb
  use mod_bas, only : znuc, pmat
  use mod_const, only : zero, iunit
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
  complex(c_double_complex), intent(out) :: d1mat(1:(2*ndvr+1), 1:(nrad-1))
  integer(c_int) :: ifun, irad, jrad, jll, jul, ij, jb1, ctdep1, ctdep2

  d1mat(1:(2*ndvr+1), 1:(nrad-1)) = zero

  do irad = 1, nrad - 1
     jll = max(1,        irad - ndvr)
     jul = min(nrad - 1, irad + ndvr)
     do jrad = jll, jul
        jb1 = ndvr + 1 + jrad - irad
        d1mat(jb1, irad) = pmat(jb1, irad)
     end do
  end do


end subroutine bas_gen_d1mat
!######################################################################


