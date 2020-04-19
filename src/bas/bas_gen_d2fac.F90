!######################################################################
subroutine bas_gen_d2fac()

  use, intrinsic :: iso_c_binding
  use mod_sph, only : mmax2, lmax2
  use mod_rad, only : nrad, xrad, wrad, ecs_flag, cxrad, inf_range, exp_factor, irad_inf
  use mod_control, only : fedvr_normalized, exact3j
  use mod_const, only : one, half, four, PI
  use mod_bas, only : bas_d2fac1, bas_d2fac2, bas_d2invr, bas_d2rpl0, bas_d2rpl1, bas_d2crpl1

  implicit none
  integer(c_long) :: l, m, irad, dim
  real(c_double) :: r2pi
! Orimo_ECS
  real(c_double) :: tmp_exp_factor
! Orimo_ECS

  dim = nrad - 1

! Orimo_ECS
  if (inf_range) then
     tmp_exp_factor = exp_factor
  else
     tmp_exp_factor = 0
  end if
! Orimo_ECS

  if (.not.exact3j) then
     ! factor r2pi comes from separately normalized l and m functions
     r2pi = half / PI
  else
!NOTE HERE. SEE ALSO hprod_mkrho2_x3j and hprod_mfprodx3j
     r2pi = 1d0
!     r2pi = half / PI
!NOTE HERE. SEE ALSO hprod_mkrho2_x3j and hprod_mfprodx3j
  end if

  if (fedvr_normalized) then
     do l = 0, lmax2
!        bas_d2fac1(l) = one / xrad(dim + 1) ** (2 * l + 1)
        do irad = 1, dim
           bas_d2rpl0(irad, l) = xrad(irad) ** l * r2pi
           bas_d2rpl1(irad, l) = xrad(irad) ** (l + 1) * sqrt(wrad(irad))
           bas_d2invr(irad, l) = (2 * l + 1) / (xrad(irad) * sqrt(wrad(irad))) * r2pi
           bas_d2fac2(irad, l) = four * pi / (2 * l + 1) / (xrad(irad) * sqrt(wrad(irad)))
        end do
! Orimo_ECS
        if (inf_range) then
           do irad = irad_inf, dim
              bas_d2invr(irad, l) = & 
                   & (2 * l + 1) & 
                   & / (xrad(irad) * sqrt(wrad(irad)) * exp((xrad(irad)-xrad(irad_inf)) * tmp_exp_factor)) & 
                   & * r2pi
           end do
        end if
        if (ecs_flag == 1) then
           do irad = 1, dim
              bas_d2crpl1(irad, l) = four * pi / (2 * l + 1) / cxrad(irad) ** (l+1)
           end do
        end if
! Orimo_ECS
     end do
  else
! Orimo_ECS
     if (ecs_flag == 1) then
        stop 'bas_gen_d2fac : not yet implemented for ECS with unnormalized fedvr'
     end if
! Orimo_ECS
     do l = 0, lmax2
!        bas_d2fac1(l) = one / xrad(dim + 1) ** (2 * l + 1)
        do irad = 1, dim
           bas_d2rpl0(irad, l) = xrad(irad) ** l * r2pi * wrad(irad)
           bas_d2rpl1(irad, l) = xrad(irad) ** (l + 1)
           bas_d2invr(irad, l) = (2 * l + 1) / xrad(irad) * wrad(irad) * r2pi
           bas_d2fac2(irad, l) = four * pi / (2 * l + 1) / xrad(irad) * wrad(irad)
        end do
     end do
  end if

end subroutine bas_gen_d2fac
!######################################################################
