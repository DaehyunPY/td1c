!######################################################################
subroutine hprod_mkrhoxc(rho1, rhoxc)
!
! rho1:  active density in the spherical harmonics basis
! rhoxc: total density on the angular grid
!
  use, intrinsic :: iso_c_binding
  use mod_hprod, only : orb, orbg
  use mod_ormas, only : nfcore, nocc, nelact, neltot
  use mod_sph, only : lmax1, lmax2, nlat, wlat
  use mod_rad, only : nrad, xrad, wrad
  use mod_const, only : zero, two, czero, pi

  implicit none
  complex(c_double_complex), intent(out) :: rho1(1:(nrad-1), 0:lmax2)
  real(c_double), intent(out) :: rhoxc(1:(nrad-1), 1:nlat)

  real(c_double) :: nume
  integer(c_int) :: llr, ulr, ifun, l, ilat, irad
  complex(c_double_complex), allocatable :: rhog(:,:)
  complex(c_double_complex), allocatable :: rhos(:,:)

  allocate(rhog(1:(nrad-1), 1:nlat))
  allocate(rhos(1:(nrad-1), 0:lmax2))
  rhog(1:(nrad-1), 1:nlat) = czero
  rhoxc(1:(nrad-1), 1:nlat) = zero

  !$omp parallel default(shared) private(llr, ulr)
  !###########################
  call util_omp_disp(1, nrad - 1, llr, ulr)
  do ifun = nfcore + 1, nocc
     do ilat = 1, nlat
        do irad = llr, ulr
           rhog(irad, ilat) = rhog(irad, ilat) + conjg(orbg(irad, ilat, ifun)) &
                                                     * orbg(irad, ilat, ifun)
        end do
     end do
  end do
  do ifun = 1, nocc
     do ilat = 1, nlat
        do irad = llr, ulr
           rhoxc(irad, ilat) = rhoxc(irad, ilat) + dble(conjg(orbg(irad, ilat, ifun)) &
                                                            * orbg(irad, ilat, ifun))
        end do
     end do
  end do

  do ilat = 1, nlat
     do irad = llr, ulr
        rhog(irad, ilat) = two * rhog(irad, ilat)
!       rhoxc(irad, ilat) =  two * rhoxc(irad, ilat) / (xrad(irad) ** two * wrad(irad))
        rhoxc(irad, ilat) =  two * rhoxc(irad, ilat) / (xrad(irad) ** two * wrad(irad)) / (two * pi)
!       rhoxc(irad, ilat) =  two * rhoxc(irad, ilat) / (xrad(irad) ** two) / (two * pi)
     end do
  end do
  !###########################
  !$omp end parallel

!debug
!  nume = zero
!  do ilat = 1, nlat
!     do irad = 1, nrad - 1
!        nume = nume + dble(rhog(irad, ilat)) * wlat(ilat)
!!        write(6, "(2i5,2f20.10)") ilat, irad, rhog(irad, ilat)
!     end do
!  end do
!  write(6, "('hprod_mkrhoxc: number of electrons from rhog  = ', f20.10, i10)") nume, neltot(3)
!  nume = zero
!  do ilat = 1, nlat
!     do irad = 1, nrad - 1
!        nume = nume + rhoxc(irad, ilat) * wlat(ilat) * xrad(irad) * xrad(irad) * wrad(irad) * (two * pi)
!     end do
!  end do
!  write(6, "('hprod_mkrhoxc: number of electrons from rhoxc = ', f20.10, i10)") nume, neltot(3)
!debug

  call bas_ang2sph2one(0, rhog, rhos)
  rho1(1:(nrad-1), 0:lmax2) = rhos(1:(nrad-1), 0:lmax2)

!debug
!  nume = zero
!  do l = 0, lmax2
!     do irad = 1, nrad - 1
!        nume = nume + dble(rho1(irad, l)) * sqrt(two)
!     end do
!  end do
!  write(6, "('hprod_mkrhoxc: number of electrons from rho1  = ', f20.10, i10)") nume, nelact(3)

!  nume = zero
!  do l = 0, lmax1
!     do irad = 1, nrad - 1
!        nume = nume + dble(conjg(orb(irad, l, 1)) * orb(irad, l, 1)) * two
!     end do
!  end do
!  write(6, "('hprod_mkrhoxc: number of electrons from orb   = ', f20.10, i10)") nume, neltot(3)
!
!  nume = zero
!  do ilat = 1, nlat
!     do irad = 1, nrad - 1
!        nume = nume + dble(conjg(orbg(irad, ilat, 1)) * orbg(irad, ilat, 1)) * wlat(ilat) * two
!     end do
!  end do
!  write(6, "('hprod_mkrhoxc: number of electrons from orbg  = ', f20.10, i10)") nume, neltot(3)
!
!  call bas_sph2ang2one(0, rhos, rhog)
!  nume = zero
!  do ilat = 1, nlat
!     do irad = 1, nrad - 1
!        nume = nume + dble(rhog(irad, ilat)) * wlat(ilat)
!!        write(6, "(2i5,2f20.10)") ilat, irad, rhog(irad, ilat)
!     end do
!  end do
!  write(6, "('hprod_mkrhoxc: number of electrons from rhog  = ', f20.10, i10)") nume, neltot(3)
!debug

  deallocate(rhos)
  deallocate(rhog)

end subroutine hprod_mkrhoxc
!######################################################################
