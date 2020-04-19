!######################################################################
subroutine bas_gen_pzfac()

  use, intrinsic :: iso_c_binding
  use mod_const, only : one, iunit
  use mod_sph, only : lmax1, mmax1
  use mod_control, only : fedvr_normalized
  use mod_bas, only : alph_lm, bas_pzfac1, bas_pzfac2
  use mod_rad, only : nrad, xrad, wrad, cxrad, ecs_flag

  implicit none
  complex(c_double_complex) :: tmp
  integer(c_int) :: l, m, irad

  if (ecs_flag == 0) then
     if (fedvr_normalized) then
        do m = -mmax1, mmax1
           do l = abs(m), lmax1
              bas_pzfac1(l, m) = -iunit * alph_lm(l, m)
              tmp = -iunit * alph_lm(l, m) * (dble(l) + one)
              do irad = 1, nrad - 1
                 bas_pzfac2(irad, l, m) = tmp / xrad(irad)
              end do
           end do
        end do
     else
        do m = -mmax1, mmax1
           do l = abs(m), lmax1
              bas_pzfac1(l, m) = -iunit * alph_lm(l, m)
              tmp = -iunit * alph_lm(l, m) * (dble(l) + one)
              do irad = 1, nrad - 1
                 bas_pzfac2(irad, l, m) = tmp / xrad(irad) * wrad(irad)
              end do
           end do
        end do
     end if
  else if(ecs_flag == 1) then
     if (fedvr_normalized) then
        do m = -mmax1, mmax1
           do l = abs(m), lmax1
              bas_pzfac1(l, m) = -iunit * alph_lm(l, m)
              tmp = -iunit * alph_lm(l, m) * (dble(l) + one)
              do irad = 1, nrad - 1
                 bas_pzfac2(irad, l, m) = tmp / cxrad(irad)
              end do
           end do
        end do
     else
        write(*,*) "ALERT : ECS with not fedvr_normalized possibly does not work."
        do m = -mmax1, mmax1
           do l = abs(m), lmax1
              bas_pzfac1(l, m) = -iunit * alph_lm(l, m)
              tmp = -iunit * alph_lm(l, m) * (dble(l) + one)
              do irad = 1, nrad - 1
                 !!=========================================================
                 !! to use ECS here, it is required to check what wights should be used here.
                 bas_pzfac2(irad, l, m) = tmp / cxrad(irad) * wrad(irad)
                 !!=========================================================
              end do
           end do
        end do
     end if
  end if
  
end subroutine bas_gen_pzfac
!######################################################################
