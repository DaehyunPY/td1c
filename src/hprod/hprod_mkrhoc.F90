!######################################################################
subroutine hprod_mkrhoc()

  use, intrinsic :: iso_c_binding
  use mod_sph, only : nlat
  use mod_hprod, only : orbg, rho1
  use mod_ormas, only : nfcore, ncore
  use mod_rad, only : nrad, xrad, wrad
  use mod_const, only : czero, ctwo, two, pi

  implicit none
  integer(c_int) :: llr, ulr, ifun, ilat, irad

  rho1(1:(nrad-1), 1:nlat) = czero
  if (ncore == 0) return

  !$omp parallel default(shared) private(llr, ulr)
  !###########################
  call util_omp_disp(1, nrad - 1, llr, ulr)
  do ifun = nfcore + 1, ncore
     do ilat = 1, nlat
        do irad = llr, ulr
           rho1(irad, ilat) = rho1(irad, ilat) + conjg(orbg(irad, ilat, ifun)) &
                                                     * orbg(irad, ilat, ifun)
        end do
     end do
  end do

  do ilat = 1, nlat
     do irad = llr, ulr
!       rho1(irad, ilat) =  ctwo * rho1(irad, ilat) / (xrad(irad) * xrad(irad) * wrad(irad))
!       rho1(irad, ilat) =  ctwo * rho1(irad, ilat) / (xrad(irad) * xrad(irad) * wrad(irad)) / (two * pi)
        rho1(irad, ilat) =  ctwo * rho1(irad, ilat) / (xrad(irad) * xrad(irad)) / (two * pi)
     end do
  end do
  !###########################
  !$omp end parallel

!debug  do irad = 1, nrad - 1
!debug     write(6, "(i5, 5f20.10)") irad, xrad(irad), &
!debug          dble(orbg(irad,1,1) / sqrt(wrad(irad))), &
!debug          dble(orbg(irad,1,2) / sqrt(wrad(irad))), &
!debug          dble(rho1(irad, 1)), &
!debug          dble(-(3.d+0/3.14d+0) ** (1.d+0/3.d+0) * rho1(irad, 1) ** (1.d+0/3.d+0))
!debug  end do
!debug  stop

end subroutine hprod_mkrhoc
!######################################################################
