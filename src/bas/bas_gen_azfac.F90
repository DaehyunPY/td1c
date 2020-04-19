!######################################################################
subroutine bas_gen_azfac()

  use, intrinsic :: iso_c_binding
  use mod_const, only : two, pi
  use mod_sph, only : lmax1, mmax1
  use mod_rad, only : nrad, xrad, wrad
  use mod_bas, only : znuc, alph_lm, bas_azfac
  use mod_control, only : fedvr_normalized, SAE, PSP, psp_type
  use mod_bas, only : pp_znuc, pp_rloc, pp_cloc

  implicit none
  complex(c_double_complex) :: aazfac, tmp
  integer(c_long) :: l, m, irad

  real(c_double), parameter ::  a1 = +1.231d+0
  real(c_double), parameter ::  a2 = +0.662d+0
  real(c_double), parameter ::  a3 = -1.325d+0
  real(c_double), parameter ::  a4 = +1.236d+0
  real(c_double), parameter ::  a5 = -0.231d+0
  real(c_double), parameter ::  a6 = +0.480d+0

  bas_azfac = 0d0

  if (fedvr_normalized) then
     do m = -mmax1, mmax1
        do l = abs(m), lmax1
           if (PSP .and. psp_type==1) then
              aazfac = -alph_lm(l, m)
              do irad = 1, nrad - 1
                 bas_azfac(irad, l, m) = aazfac &
                    * (pp_znuc/(xrad(irad))**2     *erf(  xrad(irad)/(sqrt(2d0)*pp_rloc)    ) &
                     - pp_znuc/(xrad(irad)*pp_rloc)*exp(-(xrad(irad)/(sqrt(2d0)*pp_rloc))**2)*sqrt(2.0/pi) &
                     - xrad(irad)/pp_rloc**2*exp(-5d-1*(xrad(irad)/pp_rloc)**2) &
                       *(pp_cloc(1) &
                       + pp_cloc(2)*(xrad(irad)/pp_rloc)**2 &
                       + pp_cloc(3)*(xrad(irad)/pp_rloc)**4 &
                       + pp_cloc(4)*(xrad(irad)/pp_rloc)**6) &
                     + exp(-5d-1*(xrad(irad)/pp_rloc)**2)/pp_rloc &
                       *(2d0*pp_cloc(2)*(xrad(irad)/pp_rloc) &
                       + 4d0*pp_cloc(3)*(xrad(irad)/pp_rloc)**3 &
                       + 6d0*pp_cloc(4)*(xrad(irad)/pp_rloc)**5))
              end do
           else if (SAE) then
              aazfac = -alph_lm(l, m)
              do irad = 1, nrad - 1
                 bas_azfac(irad, l, m) = aazfac / (xrad(irad))**two &
                    * ( 1.D+0 &
                      + a1 * exp(-a2*xrad(irad)) &
                      + a3 * exp(-a4*xrad(irad)) * xrad(irad) &
                      + a5 * exp(-a6*xrad(irad)) &
                      + a1 * a2 * exp(-a2*xrad(irad)) &
                      + a3 * a4 * exp(-a4*xrad(irad)) * xrad(irad) &
                      + a5 * a6 * exp(-a6*xrad(irad)) &
                      - a3 * exp(-a4*xrad(irad)))
              end do
           else
              aazfac = -alph_lm(l, m) * znuc
              do irad = 1, nrad - 1
                 bas_azfac(irad, l, m) = aazfac / (xrad(irad))**two
              end do
           end if
        end do
     end do
  else
     do m = -mmax1, mmax1
        do l = abs(m), lmax1
           aazfac = -alph_lm(l, m) * znuc
           do irad = 1, nrad - 1
              bas_azfac(irad, l, m) = aazfac / (xrad(irad))**two * wrad(irad)
           end do
        end do
     end do
  end if

end subroutine bas_gen_azfac
!######################################################################
