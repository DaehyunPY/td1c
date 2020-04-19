!######################################################################
real(c_double) function hprod_ene_fcx()

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nradfc
  use mod_ormas, only : nfcore2
  use mod_hprod, only : rhofc, v2xfc
  use mod_const, only : zero, two, quart

  implicit none
  integer(c_long), external :: util_omp_nproc
  integer(c_long), external :: util_omp_iproc
  integer(c_long) :: iproc, nproc, irad, ilat, llr, ulr
  real(c_double) :: etmp1, enep

  if (nfcore2 == 0) then
     hprod_ene_fcx = zero
     return
  end if

  enep = zero
  nproc = util_omp_nproc()
  !$omp parallel default(shared) private(iproc, etmp1, llr, ulr) reduction(+:enep)
  !###########################
  iproc = util_omp_iproc()
  call util_omp_disp(1, nradfc, llr, ulr)
  etmp1 = zero
  do irad = llr, ulr
     etmp1 = etmp1 + rhofc(irad) * v2xfc(irad)
  end do
  enep = enep + etmp1
  !###########################
  !$omp end parallel

  hprod_ene_fcx = two * enep * quart
!DEBUG
!  write(6, "('ene_fcore_2 = ', f20.10)") hprod_ene_fcx
!DEBUG

end function hprod_ene_fcx
!######################################################################
