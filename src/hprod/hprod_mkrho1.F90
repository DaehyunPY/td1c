!######################################################################
subroutine hprod_mkrho1(rho1)

  use, intrinsic :: iso_c_binding
  use mod_hprod, only : orbg
  use mod_ormas, only : nocc, neltot
  use mod_sph, only : nlat, wlat
  use mod_rad, only : nrad, xrad, wrad
  use mod_const, only : zero, two, pi

  implicit none
  real(c_double), intent(out) :: rho1(1:(nrad-1), 1:nlat)
  real(c_double) :: nume
  integer(c_int) :: llr, ulr, ifun, ilat, irad

  rho1(1:(nrad-1), 1:nlat) = zero

  !$omp parallel default(shared) private(llr, ulr)
  !###########################
  call util_omp_disp(1, nrad - 1, llr, ulr)
  do ifun = 1, nocc
     do ilat = 1, nlat
        do irad = llr, ulr
           rho1(irad, ilat) = rho1(irad, ilat) + dble(conjg(orbg(irad, ilat, ifun)) &
                                                          * orbg(irad, ilat, ifun))
        end do
     end do
  end do

  do ilat = 1, nlat
     do irad = llr, ulr
!       rho1(irad, ilat) =  two * rho1(irad, ilat) / (xrad(irad) * xrad(irad) * wrad(irad))
        rho1(irad, ilat) =  two * rho1(irad, ilat) / (xrad(irad) * xrad(irad) * wrad(irad)) / (two * pi)
!       rho1(irad, ilat) =  two * rho1(irad, ilat) / (xrad(irad) * xrad(irad)) / (two * pi)
     end do
  end do
  !###########################
  !$omp end parallel

!debug
  nume = zero
  do ilat = 1, nlat
     do irad = 1, nrad - 1
        nume = nume + rho1(irad, ilat) * wlat(ilat) * xrad(irad) * xrad(irad) * wrad(irad) * (two * pi)
        write(6, "(2i5,f20.10)") ilat, irad, rho1(irad, ilat)
     end do
  end do
  write(6, "('hprod_mkrho1: number of electrons = ', f20.10, i10)") nume, neltot(3)
!debug

!debug  do irad = 1, nrad - 1
!debug     write(6, "(i5, 5f20.10)") irad, xrad(irad), &
!debug          dble(orbg(irad,1,1) / sqrt(wrad(irad))), &
!debug          dble(orbg(irad,1,2) / sqrt(wrad(irad))), &
!debug          dble(rho1(irad, 1)), &
!debug          dble(-(3.d+0/3.14d+0) ** (1.d+0/3.d+0) * rho1(irad, 1) ** (1.d+0/3.d+0))
!debug  end do
!debug  stop

end subroutine hprod_mkrho1
!######################################################################
