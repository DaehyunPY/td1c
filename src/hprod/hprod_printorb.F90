!######################################################################
subroutine hprod_printorb(wfn)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad, xrad, wrad
  use mod_sph, only : lmax1
  use mod_bas, only : mval
  use mod_ormas, only : nfun

  implicit none
  complex(c_double_complex), intent(in) :: wfn(1:(nrad-1), 0:lmax1, 1:*)
  integer(c_int) :: ifun, irad, l

  do ifun = 1, nfun
     write(6,"('# ifun = ',i5)") ifun
     do irad = 1, nrad-1
        write(6,"(f10.5)",advance='no') xrad(irad)
        do l = 0, lmax1
           write(6,"(2f20.10)",advance='no') wfn(irad,l,ifun)/sqrt(wrad(irad))
        end do
        write(6,*)
     end do
     write(6,*)
     write(6,*)
  end do

end subroutine hprod_printorb
!######################################################################
subroutine hprod_printrho(rho)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad, xrad
  use mod_sph, only : lmax2
  use mod_bas, only : mval
  use mod_ormas, only : nfun

  implicit none
  complex(c_double_complex), intent(in) :: rho(1:(nrad-1), 0:lmax2, 1:nfun, 1:nfun)
  integer(c_int) :: ifun, jfun, irad, l

  do ifun = 1, nfun
  do jfun = 1, nfun
     write(6,"('# ifun = ',2i5)") ifun, jfun
     do irad = 1, nrad-1
        write(6,"(f10.5)",advance='no') xrad(irad)
        do l = 0, lmax2
           write(6,"(2f20.10)",advance='no') rho(irad,l,ifun,jfun)
        end do
        write(6,*)
     end do
     write(6,*)
     write(6,*)
  end do
  end do

end subroutine hprod_printrho
!######################################################################
subroutine hprod_printrhog(rhog)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad, xrad
  use mod_sph, only : nlat
  use mod_bas, only : mval
  use mod_ormas, only : nfun

  implicit none
  complex(c_double_complex), intent(in) :: rhog(1:(nrad-1), 1:nlat, 1:nfun, 1:nfun)
  integer(c_int) :: ifun, jfun, irad, ilat

  do ifun = 1, nfun
  do jfun = 1, nfun
     write(6,"('# ifun = ',2i5)") ifun, jfun
     do irad = 1, nrad-1
        write(6,"(f10.5)",advance='no') xrad(irad)
        do ilat = 1, nlat
           write(6,"(2f20.10)",advance='no') rhog(irad,ilat,ifun,jfun)
        end do
        write(6,*)
     end do
     write(6,*)
     write(6,*)
  end do
  end do

end subroutine hprod_printrhog
!######################################################################
