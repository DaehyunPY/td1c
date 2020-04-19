!######################################################################
subroutine calc_orb(sirad, orb, norb, dorb)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad
  use mod_bas, only : nbas
  use mod_ormas, only : nfun , nstr_alph, nstr_beta
  use mod_sph, only : lmax1

  implicit none
  integer(c_long), intent(in) :: sirad
  complex(c_double_complex), intent(in) :: orb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(out) :: norb(0:lmax1, 1:nfun)
  complex(c_double_complex), intent(out) :: dorb(0:lmax1, 1:nfun)

  integer(c_long) :: lll, ull
  
  !$omp parallel default(shared) private(lll, ull)
  !###########################
  call util_omp_disp(0, lmax1, lll, ull)
  call calc_orbp(lll, ull, sirad, orb, norb, dorb)
  !###########################
  !$omp end parallel

  !DEBUG
  !write(6, "('# calc_orb_: done', f20.10)") nfun
  !DEBUG



end subroutine calc_orb
!######################################################################
subroutine calc_orbp(lll, ull, sirad, orb, norb, dorb)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nfun, nfcore, ncore
  use mod_rad, only : nrad, ndvr, xrad, cxrad, wrad, cwrad, irad_ecs, mapf, mapb
  use mod_sph, only : lmax1
  use mod_const, only : czero, ctwo, iunit
  use mod_bas, only : d1mat, d2mat

  implicit none
  integer(c_long), intent(in) :: lll, ull
  integer(c_long), intent(in) :: sirad
  complex(c_double_complex), intent(in) :: orb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(out) :: norb(0:lmax1, 1:nfun)
  complex(c_double_complex), intent(out) :: dorb(0:lmax1, 1:nfun)

  integer(c_long) :: ifun, l, irad, jrad, jb1, kb1
  integer(c_long) :: jb, ife, iloc, jloc

  complex(c_double_complex) :: d1tmp


  do ifun = nfcore + 1, nfun

     do l = lll, ull     
        do irad = sirad, sirad

           ife = mapb(irad)
           iloc = irad - mapf(ife);           
           
           d1tmp = czero
           ! non bridge elements
           do jloc = 0, ndvr
              jrad = mapf(ife) + jloc;
              if (jrad > 0 .and. jrad < nrad) then
                 jb = ndvr + 1 + jrad - irad
                 d1tmp = d1tmp + orb(jrad, l, ifun) * d1mat(jb, irad)
              end if
           end do

           ! bridge functions
           if (iloc == 0 .and. ife /= 0) then
              do jloc = 0, ndvr - 1
                 jrad = mapf(ife - 1) + jloc;
                 if (jrad > 0 .and. jrad < nrad) then
                    jb = ndvr + 1 + jrad - irad
                    d1tmp = d1tmp + orb(jrad, l, ifun) * d1mat(jb, irad)
                 end if
              end do
           end if

           norb(l, ifun) = orb(sirad, l , ifun) / (cxrad(sirad) * sqrt(cwrad(sirad)))
           dorb(l, ifun) = (d1tmp - orb(sirad, l, ifun) / cxrad(sirad)) / (cxrad(sirad) * sqrt(cwrad(sirad)))

        end do
     end do
  end do


end subroutine calc_orbp
!######################################################################
