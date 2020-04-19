!######################################################################
! test FC use nfcore_tdcis instead of nfcore
! test FC use nradgs instead of nrad
!######################################################################
subroutine hprod_mkrho2_tdcis()

  use, intrinsic :: iso_c_binding
  use mod_hprod, only : orb, orbg, rho2, v2ang
  use mod_hprod, only : orb0rotg, orbg
  use mod_rad, only : nrad, nradfc
  use mod_sph, only : nlat
  use mod_ormas, only : nfun, nfcore
  use mod_ormas, only : nfcore_tdcis
  use mod_rad, only : nradgs
  use mod_rad, only : ecs_flag, irad_ecs

  implicit none
  integer(c_int) :: llr, ulr, ifun, jfun, ilat, irad, nradmax

  if (ecs_flag == 0) then
     nradmax = nrad - 1
  else 
     nradmax = irad_ecs - 1
  end if

  !$omp parallel default(shared) private(llr, ulr)
  !###########################
  call util_omp_disp(1, irad_ecs, llr, ulr) 

  ! does not work well 
  ! call util_omp_disp(1, nradgs - 1, llr, ulr)
  ! does not work well 

  ! do ifun = nfcore + 1, nfun
  !    do jfun = nfcore + 1, nfun
  do ifun = nfcore_tdcis + 1, nfun
     do jfun = nfcore_tdcis + 1, nfun
        do ilat = 1, nlat
           do irad = llr, ulr
             v2ang(irad, ilat, jfun, ifun) = conjg(orb0rotg(irad, ilat, jfun)) &
                                             * orbg(irad, ilat, ifun)
           end do
        end do
     end do
  end do
  !###########################
  !$omp end parallel
  call bas_ang2sph2(v2ang, rho2);

end subroutine hprod_mkrho2_tdcis
!######################################################################
! this subroutine is called only from hprod_set_wfn0
subroutine hprod_mkrho2_tdcis_init()

  use, intrinsic :: iso_c_binding
  use mod_hprod, only : orb, orbg, rho2, v2ang
  use mod_hprod, only : orb0rotg, orbg
  use mod_rad, only : nrad, nradfc
  use mod_sph, only : nlat
  use mod_ormas, only : nfun
  use mod_rad, only : nradgs
  use mod_rad, only : ecs_flag, irad_ecs

  implicit none
  integer(c_int) :: llr, ulr, ifun, jfun, ilat, irad, nradmax

  if (ecs_flag == 0) then
     nradmax = nrad - 1
  else 
     nradmax = irad_ecs - 1
  end if

  !$omp parallel default(shared) private(llr, ulr)
  !###########################
  call util_omp_disp(1, nradmax, llr, ulr)

  ! does not work well 
  ! call util_omp_disp(1, nradgs - 1, llr, ulr)
  ! does not work well 

  do ifun = 1, nfun
     do jfun = 1, nfun
        do ilat = 1, nlat
           do irad = llr, ulr
             v2ang(irad, ilat, jfun, ifun) = conjg(orb0rotg(irad, ilat, jfun)) &
                                             * orbg(irad, ilat, ifun)
           end do
        end do
     end do
  end do
  !###########################
  !$omp end parallel
  call bas_ang2sph2(v2ang, rho2);

end subroutine hprod_mkrho2_tdcis_init
!######################################################################
