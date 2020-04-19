!######################################################################
subroutine hprod_orbin_tdcis_init()

  use, intrinsic :: iso_c_binding
  use mod_bas, only : nbas
  use mod_ormas, only : nfun
  use mod_hprod, only : orb

  implicit none

  ! write(6, *) nbas
  call zclear_omp(nbas*nfun, orb)

end subroutine hprod_orbin_tdcis_init
!######################################################################
subroutine hprod_orbin_tdcis(lfield)

  use, intrinsic :: iso_c_binding
  ! debug
  use mod_rad, only : nrad
  use mod_sph, only : nlat
  ! debug
  use mod_bas, only : nbas
  use mod_control, only : igauge, tdcis_rvg
  use mod_ormas, only : nfcore, nfun
  use mod_hprod, only : orb0, orb0rot, orb0rotg

  implicit none
  real(c_double), intent(in) :: lfield(1:3, 1:3)

  ! debug
  integer(c_int) :: size
  size = (nrad - 1)*nlat*nfun
  ! debug

  if (tdcis_rvg) then
     ! rVG
     call hprod_fcorbg_tdcis(lfield)
     call bas_ang2sph1(orb0rotg, orb0rot)
     ! debug
     ! call util_print_vec(size, orb0rotg, 'orbin.orb0rotg')
     ! debug
  else
     ! LG or VG
     orb0rot = orb0
     call bas_sph2ang1(orb0rot, orb0rotg)
     ! debug
     ! call util_print_vec(size, orb0rotg, 'orbin.orb0rotg')
     ! debug
  end if

end subroutine hprod_orbin_tdcis
!######################################################################
subroutine hprod_fcorbg_tdcis(lfield)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : xrad, nrad, nradfc
  use mod_sph, only : lmax1, nlat, cost
  use mod_ormas, only : nfcore, nfun
  use mod_hprod, only : orb0, orb0rotg
  use mod_const, only : iunit

  implicit none
  real(c_double), intent(in) :: lfield(1:3, 1:3)

  integer(c_int) :: ifun
  integer(c_int) :: llr, ulr, irad, ilat
  complex(c_double_complex) :: afield, efield
  complex(c_double_complex) :: faclv

  ! debug
  integer(c_int):: size
  size = (nrad - 1)*nlat*nfun
  ! debug

  efield = lfield(3, 2) ! electric field
  afield = lfield(3, 3) ! vector potential

  call bas_sph2ang1(orb0, orb0rotg)


  ! debug
  ! call util_print_vec(size, orb0rotg, 'fcorbg.orb0rotg')
  ! debug

  !$omp parallel default(shared) private(faclv, llr, ulr)
  !###########################
  call util_omp_disp(1, nradfc, llr, ulr)
  do ilat = 1, nlat
     do irad = llr, ulr
        faclv = exp(-iunit * afield * xrad(irad) * cost(ilat) )
        ! faclv = exp(-iunit * lfield(3, 3) * xrad(irad) * cost(ilat))
        do ifun = 1, nfun
           orb0rotg(irad, ilat, ifun) = orb0rotg(irad, ilat, ifun) * faclv
        end do
     end do
  end do
  !###########################
  !$omp end parallel

end subroutine hprod_fcorbg_tdcis
!######################################################################
