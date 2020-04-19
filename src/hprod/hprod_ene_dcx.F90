!######################################################################
real(c_double) function hprod_ene_dcx()

  use, intrinsic :: iso_c_binding
  use mod_hprod, only : rho1
  use mod_sph, only : nlat, wlat
  use mod_rad, only : nrad, xrad, wrad
  use mod_ormas, only : ndcore
  use mod_const, only : zero, two, four, third, quart, fac_lda, pi

  implicit none
  integer(c_int), external :: util_omp_nproc
  integer(c_int), external :: util_omp_iproc
  integer(c_int) :: iproc, nproc, irad, ilat, llr, ulr
  real(c_double) :: etmp1, etmp2, enep

  if (ndcore == 0) then
     hprod_ene_dcx = zero
     return
  end if

  enep = zero
  nproc = util_omp_nproc()
  !$omp parallel default(shared) private(iproc, etmp1, etmp2, llr, ulr) reduction(+:enep)
  !###########################
  iproc = util_omp_iproc()
  call util_omp_disp(1, nrad - 1, llr, ulr)
  etmp1 = zero
  do ilat = 1, nlat
     etmp2 = zero
     do irad = llr, ulr
        etmp2 = etmp2 + wrad(irad) * xrad(irad)**two * rho1(irad, ilat)**(four*third)
     end do
     etmp1 = etmp1 + etmp2 * wlat(ilat)
  end do
  enep = enep + etmp1
  !###########################
  !$omp end parallel

! hprod_ene_dcx = zero
! hprod_ene_dcx = -enep * fac_lda * quart
  hprod_ene_dcx = - enep * fac_lda * quart * two * pi
! hprod_ene_dcx = +enep * fac_lda * quart
! hprod_ene_dcx = -enep * fac_lda * 1.d+0 / 8.d+0
! hprod_ene_dcx = +enep * fac_lda * 1.d+0 / 8.d+0
!DEBUG
!  write(6, "('ene_dcore_2 = ', f20.10)") hprod_ene_dcx
!DEBUG

end function hprod_ene_dcx
!######################################################################
