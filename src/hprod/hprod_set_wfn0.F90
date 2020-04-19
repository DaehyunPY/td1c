!#######################################################################
subroutine hprod_set_wfn0(wfn, cic, v2xfc_type)

  use, intrinsic :: iso_c_binding
  use mod_control, only : name
  use mod_bas, only : nbas
  use mod_rad, only : nrad
  use mod_const, only : zero, czero
  use mod_hprod, only : orb0, cic0, v2xfc, projhigh
  use mod_ormas, only : nfun, lcic, nfcore, nfcore2
! tdcis-teramura
  use mod_control, only : istdcis
  use mod_rad, only : xrad
  use mod_ormas, only : thrdet
  use mod_hprod, only : ridm_tdcis
  use mod_hprod, only : tdcis_eig, v2ang0, orb0rot, orb0rotg, h0orb0, orb, orbg, v2ang
! tdcis-teramura

  implicit none
  complex(c_double_complex), intent(in) :: wfn(1:*)
  complex(c_double_complex), intent(in) :: cic(1:*)
  character(*), intent(in) :: v2xfc_type

  character(len=256) :: fname
  integer(c_int) :: irad, nrad0

! tdcis-teramura
  integer(c_int) :: ifun
  real(c_double) :: lfield0(1:3,1:3), R0
  complex(c_double_complex) :: ovlpin(1:nfun, 1:nfun)
  complex(c_double_complex), external :: util_det
! tdcis-teramura

  call zcopy_omp(nbas*nfun, wfn, orb0)
  call zcopy_omp(lcic, cic, cic0)

! tdcis-teramura
  if (istdcis) then
     ! GS energy
     lfield0 = 0d0
     call hprod_orbene(lfield0, wfn, cic, tdcis_eig)
     write(6, "('hprod_set_wfn0: tdcis_eig:') ")
     write(6, "(f20.10)") tdcis_eig(1:nfun)

!debug
!     stop 'here 1'
!debug

     ! locality
     R0 = 20d0
     call hprod_mkovlp(R0, orb0, orb0, ovlpin)
     write(6, "('hprod_set_wfn0: det_S0 = ', 2es20.10)") util_det(nfun, thrdet , ovlpin)

     ! debug
     ridm_tdcis = czero
     call hprod_mkovlp(xrad(nrad), orb0, orb0, ridm_tdcis)
     write(6 , "('hprod_set_wfn0: debug ')") 
     do ifun = 1, nfun
        write(6 , "(4es20.10)") ridm_tdcis(ifun, ifun), ovlpin(ifun, ifun)
     end do
     ! debug

     ! orb0rot
     call hprod_orbin_tdcis(lfield0)
     orb = orb0rot
     orbg = orb0rotg

     ! h0orb0 
     call hprod_tprod_all(orb0rot, h0orb0)
     
     ! v2ang0
     call hprod_mkrho2_tdcis_init
     call hprod_mkv2mf_tdcis_init
     v2ang0 = v2ang

!debug
!     stop 'here 2'
!debug

     ! debug
     ! call util_print_vec(osize, h0orb0, 'setwfn.h0orb0')
     ! call util_print_vec(osize, orb, 'setwfn.cis.orb')
     ! call util_print_vec(osize, orb0, 'setwfn.cis.orb0')
     ! call util_print_vec(osize, orb0rot, 'setwfn.cis.orb0rot')
     ! call util_print_vec(v2size, v2ang0, 'setwfn.cis.v2ang')
     ! stop
     ! debug
  end if
  ! call util_print_vec(osize, orb0, 'setwfn.orb0')
  ! stop
! tdcis-teramura

!NEW
  if (projhigh) then
     call hprod_getprojhigh
  end if
!  stop 'for debug after hprod_getprojhigh'
!NEW

  if (nfcore2 > 0 .and. trim(v2xfc_type) == "read") then
     fname = trim(name)//".v2xfc"
     open(unit=1,file=trim(fname),status='old',form='formatted')
     read(1, *) nrad0
!c++ do irad = 1, nrad0 - 1
     do irad = 1, nrad0
        read(1, *) v2xfc(irad)
     end do
!c++ v2xfc(nrad0:nrad-1) = czero
     v2xfc((nrad0+1):(nrad-1)) = czero
     close(unit=1)
  end if

  call hprod_v2fc(wfn, cic, v2xfc_type);

  if (nfcore2 > 0) then
     fname = trim(name)//".v2xfc"
     open(unit=1,file=trim(fname),status='unknown',form='formatted')
     write(1, "(i10)") nrad
!c++ do irad = 1, nrad - 1
     do irad = 1, nrad
        write(1, "(2e25.15)") v2xfc(irad)
     end do
     close(unit=1)
  end if

!debug
!  stop 'here 3'
!debug

end subroutine hprod_set_wfn0
!#######################################################################
subroutine hprod_v2fc(wfn, cic, v2xfc_type)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad, ecs_flag
  use mod_bas, only : nbas
  use mod_sph, only : lmax2, nang
  use mod_const, only : zero, czero
  use mod_control, only : xfc_implicit, exact3j
  use mod_ormas, only : neltot, nfcore1, nfcore2, nfcore
  use mod_hprod, only : orb, orbg, rhofc, v2jfc, v2xfc, rho2, v2sph, v2ang, ene_fcore

  implicit none
  character(*), intent(in) :: v2xfc_type
  complex(c_double_complex), intent(in) :: wfn(*)
  complex(c_double_complex), intent(in) :: cic(*)
  real(c_double) :: lfield(1:9), enefc

  if (nfcore == 0 .or. neltot(3) < 2) then
     ene_fcore = czero
     rhofc(1:(nrad-1)) = czero ! Density of FC2+FC1
     v2jfc(1:(nrad-1)) = czero ! Hartree of FC2+FC1
     v2xfc(1:(nrad-1)) = czero ! Local-X of FC2
     return
  end if

  rhofc(1:(nrad-1)) = czero ! Density of FC2+FC1
  v2jfc(1:(nrad-1)) = czero ! Hartree of FC2+FC1
!bug_for_read v2xfc(1:(nrad-1)) = czero ! Local-X of FC2
  rho2(1:(nrad-1), 0:lmax2, 1:nfcore, 1:nfcore) = czero
  v2sph(1:(nrad-1), 0:lmax2, 1:nfcore, 1:nfcore) = czero
  v2ang(1:(nrad-1), 1:nang, 1:nfcore, 1:nfcore) = czero

  lfield(1:9) = zero
  call hprod_orbin(lfield, wfn, orb, orbg)

  if (.not.exact3j) then
     if (ecs_flag == 0) then
        call hprod_mkrhofc(orb, rhofc)   ! FC-density
        call hprod_mkv2jfc(rhofc, v2jfc) ! FC-Hartree
        call hprod_mkrho2p('fcfc', orbg, v2ang)        ! FC-FC pair density, ang
        call bas_ang2sph2_fcfc(v2ang, rho2);           ! FC-FC pair density, sph
        call hprod_mkv2mf_poisson('fcfc', rho2, v2sph) ! FC-FC pair potential, sph
        call bas_sph2ang2_fcfc(v2sph, v2ang);          ! FC-FC pair potential, ang
        call hprod_mkv2mf_herm('fcfc', v2sph, v2ang)   ! Hermitialization
     else
        call hprod_mkrhofc_ecs(orb, rhofc)   ! FC-density
        call hprod_mkv2jfc_ecs(rhofc, v2jfc) ! FC-Hartree
        call hprod_mkrho2p_ecs('fcfc', orbg, v2ang)        ! FC-FC pair density, ang
        call bas_ang2sph2_fcfc(v2ang, rho2);               ! FC-FC pair density, sph
        call hprod_mkv2mf_poisson_ecs('fcfc', rho2, v2sph) ! FC-FC pair potential, sph
        call bas_sph2ang2_fcfc_ecs(v2sph, v2ang);          ! FC-FC pair potential, ang
        !not hermitian call hprod_mkv2mf_herm('fcfc', v2sph, v2ang)   ! Hermitialization
     end if
  else
     call hprod_mkrhofc(orb, rhofc)   ! FC-density
     call hprod_mkv2jfc(rhofc, v2jfc) ! FC-Hartree
     call hprod_mkrho2p_x3j('fcfc', orb, rho2)          ! FC-FC pair density, sph
     call hprod_mkv2mf_x3j_poisson('fcfc', rho2, v2sph) ! FC-FC pair potential, sph
     call hprod_mkv2mf_x3j_herm('fcfc', v2sph)          ! Hermitialization
  end if

!  if (nfcore1 > 0) then
!     call hprod_mkrho2p('fc1fc1', orbg, v2ang)        ! FC1-FC1 pair density, ang
!     call bas_ang2sph2_fc1fc1(v2ang, rho2);          ! FC1-FC1 pair density, sph
!     call hprod_mkv2mf_poisson('fc1fc1', rho2, v2sph) ! FC1-FC1 pair potential, sph
!     call bas_sph2ang2_fc1fc1(v2sph, v2ang);         ! FC1-FC1 pair potential, ang
!     call hprod_mkv2mf_herm('fc1fc1', v2sph, v2ang)   ! Hermitialization
!  end if

  if (nfcore2 > 0) then
     if (trim(v2xfc_type) /= 'read') then
        if (trim(v2xfc_type) == 'lda') then
           call hprod_v2fc_xlda(v2xfc)
        else if (trim(v2xfc_type) == 'kli') then
           call hprod_v2fc_xkli(v2xfc)
        end if
!       call hprod_v2xfc_opt1(orbg, v2xfc)      ! 1-e orb based
        call hprod_v2xfc_opt2(cic, orbg, v2xfc) ! N-e wfn based
     end if
     if (xfc_implicit) v2jfc(1:(nrad-1)) = v2jfc(1:(nrad-1)) + v2xfc(1:(nrad-1))
  end if

  ! add FC-Hartree (and Slater) to atomic hamiltonian
  call hprod_v2fc_enefc(enefc)
  call hprod_addfc_tmat()
  ene_fcore = enefc

end subroutine hprod_v2fc
!######################################################################
subroutine hprod_mkrhofc(orb, rhofc)

  use, intrinsic :: iso_c_binding
  use mod_sph, only : lmax1
  use mod_ormas, only : nfun, nfcore
  use mod_bas, only : nval, lval, mval
  use mod_const, only : czero, one, two
  use mod_rad, only : nrad, nradfc, xrad

  implicit none
  complex(c_double_complex), intent(in) :: orb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(out) :: rhofc(1:(nrad-1))
  integer(c_int) :: llr, ulr, ifun, irad
  real(c_double) :: lfac

  rhofc(1:(nrad-1)) = czero
  if (nfcore == 0) return

  !$omp parallel default(shared) private(lfac, llr, ulr)
  !###########################
  call util_omp_disp(1, nradfc, llr, ulr)
  do ifun = 1, nfcore
     if (mval(ifun) == 0) then
        lfac = two * lval(ifun) + one
        do irad = llr, ulr
           rhofc(irad) = rhofc(irad) + conjg(orb(irad, lval(ifun), ifun)) * orb(irad, lval(ifun), ifun) * lfac
        end do
     end if
  end do
  !###########################
  !$omp end parallel

end subroutine hprod_mkrhofc
!######################################################################
subroutine hprod_mkrhofc_ecs(orb, rhofc)

  use, intrinsic :: iso_c_binding
  use mod_sph, only : lmax1
  use mod_ormas, only : nfun, nfcore
  use mod_bas, only : nval, lval, mval
  use mod_const, only : czero, one, two
  use mod_rad, only : nrad, nradfc, xrad, irad_ecs

  implicit none
  complex(c_double_complex), intent(in) :: orb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(out) :: rhofc(1:(nrad-1))
  integer(c_int) :: llr, ulr, ifun, irad, dim
  real(c_double) :: lfac

  rhofc(1:(nrad-1)) = czero
  if (nfcore == 0) return

  dim = min(nradfc, irad_ecs-1) ! <<<<<< added
  !$omp parallel default(shared) private(lfac, llr, ulr)
  !###########################
  call util_omp_disp(1, dim, llr, ulr) ! <<<<<< changed
  do ifun = 1, nfcore
     if (mval(ifun) == 0) then
        lfac = two * lval(ifun) + one
        do irad = llr, ulr
           rhofc(irad) = rhofc(irad) + conjg(orb(irad, lval(ifun), ifun)) * orb(irad, lval(ifun), ifun) * lfac
        end do
     end if
  end do
  !###########################
  !$omp end parallel

end subroutine hprod_mkrhofc_ecs
!######################################################################
subroutine hprod_v2fc_xlda(v2xfc)

  use, intrinsic :: iso_c_binding
  use mod_sph, only : lmax1
  use mod_hprod, only : orb
  use mod_ormas, only : nfun, nfcore2
  use mod_bas, only : nval, lval, mval
  use mod_rad, only : nrad, nradfc, xrad
  use mod_const, only : czero, one, two, four, pi, third, fac_lda

  implicit none
  complex(c_double_complex), intent(out) :: v2xfc(1:(nrad-1))

  real(c_double) :: lfac
  complex(c_double_complex) :: rho
  integer(c_int) :: llr, ulr, ifun, irad

  v2xfc(1:(nrad-1)) = czero
  if (nfcore2 == 0) return

  !$omp parallel default(shared) private(llr, ulr, lfac, rho)
  !###########################
  call util_omp_disp(1, nradfc, llr, ulr)
  do irad = llr, ulr
     rho = czero
     do ifun = 1, nfcore2
        if (mval(ifun) == 0) then
           lfac = two * lval(ifun) + one
           rho = rho + conjg(orb(irad, lval(ifun), ifun)) * orb(irad, lval(ifun), ifun) * lfac
        end if
     end do
     rho = rho / (four * pi * xrad(irad)**2) * two
     v2xfc(irad) = -fac_lda * rho ** third
  end do
  !###########################
  !$omp end parallel

end subroutine hprod_v2fc_xlda
!######################################################################
subroutine hprod_mkv2jfc(rhofc, v2jfc)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nfcore
  use mod_bas, only : mval, d2ll
  use mod_rad, only : ndvr, nrad, xrad
  use mod_control, only : exact3j
  use mod_const, only : one, czero, zero, iunit, pi
  use mod_bas, only : bas_d2fac1, bas_d2fac2, bas_d2invr, bas_d2rpl0, bas_d2rpl1

  implicit none
  complex(c_double_complex), intent(in) :: rhofc(1:(nrad-1))
  complex(c_double_complex), intent(out) :: v2jfc(1:(nrad-1))

  real(c_double) :: tmpr, tmpi
  real(c_double), allocatable :: rrho2(:,:)
  integer(c_int) :: dim, irad, ld, info

  v2jfc(1:(nrad-1)) = czero
  if (nfcore == 0) return

  dim = nrad - 1
  ld = ndvr + 1
  allocate(rrho2(1:dim, 1:2))

  bas_d2fac1(0) = one / xrad(dim + 1)

  tmpr = zero
  tmpi = zero
  do irad = 1, dim
     rrho2(irad, 1) = dble( rhofc(irad))
     rrho2(irad, 2) = aimag(rhofc(irad))
     tmpr = tmpr + rrho2(irad, 1) * bas_d2rpl0(irad, 0)
     tmpi = tmpi + rrho2(irad, 2) * bas_d2rpl0(irad, 0)
     rrho2(irad, 1) = rrho2(irad, 1) * bas_d2invr(irad, 0)
     rrho2(irad, 2) = rrho2(irad, 2) * bas_d2invr(irad, 0)
  end do
  tmpr = tmpr * bas_d2fac1(0)
  tmpi = tmpi * bas_d2fac1(0)
  call DPBTRS('L', dim, ndvr, 2, d2ll(1,1,0), ld, rrho2, dim, info)
  if (info /= 0) then
     write(6, "('hprod_mkv2jfc: dpbtrs bad info. ', i5)") info; stop
  end if
  do irad = 1, dim
     v2jfc(irad) = (rrho2(irad, 1) + tmpr * bas_d2rpl1(irad, 0)) * bas_d2fac2(irad, 0) &
                 + (rrho2(irad, 2) + tmpi * bas_d2rpl1(irad, 0)) * bas_d2fac2(irad, 0) * iunit
  end do

!SEE addT1  if (exact3j) then
!SEE addT1     do irad = 1, dim
!SEE addT1        v2jfc(irad) = v2jfc(irad) / (2d0*pi)
!SEE addT1     end do
!SEE addT1  end if

  deallocate(rrho2)

end subroutine hprod_mkv2jfc
!######################################################################
subroutine hprod_mkv2jfc_ecs(rhofc, v2jfc)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nfcore
  use mod_bas, only : mval, d2ll
  use mod_rad, only : ndvr, nrad, xrad, irad_ecs
  use mod_const, only : one, czero, zero, iunit
  use mod_bas, only : bas_d2fac1, bas_d2fac2, bas_d2invr, bas_d2rpl0, bas_d2rpl1, bas_d2crpl1 ! <<<<<< added

  implicit none
  complex(c_double_complex), intent(in) :: rhofc(1:(nrad-1))
  complex(c_double_complex), intent(out) :: v2jfc(1:(nrad-1))

  real(c_double) :: tmpr, tmpi
  real(c_double), allocatable :: rrho2(:,:)
  integer(c_int) :: dim, irad, ld, info

! Orimo_ECS
  complex(c_double_complex) :: charge
  integer(c_int) :: dimbc
! Sato  dimbc = min(dim, irad_ecs-1)
! Orimo_ECS

  v2jfc(1:(nrad-1)) = czero
  if (nfcore == 0) return

  dim = nrad - 1
! Sato
  dimbc = min(dim, irad_ecs-1)
! Sato
  ld = ndvr + 1
  allocate(rrho2(1:dim, 1:2))

! Orimo_ECS
  bas_d2fac1(0) = one / xrad(dimbc + 1)
  !bas_d2fac1(0) = one / xrad(dim + 1)
! Orimo_ECS
  
  tmpr = zero
  tmpi = zero
  do irad = 1, dimbc ! <<<<<<< changed
     rrho2(irad, 1) = dble( rhofc(irad))
     rrho2(irad, 2) = aimag(rhofc(irad))
     tmpr = tmpr + rrho2(irad, 1) * bas_d2rpl0(irad, 0)
     tmpi = tmpi + rrho2(irad, 2) * bas_d2rpl0(irad, 0)
     rrho2(irad, 1) = rrho2(irad, 1) * bas_d2invr(irad, 0)
     rrho2(irad, 2) = rrho2(irad, 2) * bas_d2invr(irad, 0)
  end do
  tmpr = tmpr * bas_d2fac1(0)
  tmpi = tmpi * bas_d2fac1(0)
! Orimo_ECS
  call DPBTRS('L', dimbc, ndvr, 2, d2ll(1,1,0), ld, rrho2, dim, info)
! Orimo_ECS
  if (info /= 0) then
     write(6, "('hprod_mkv2jfc: dpbtrs bad info. ', i5)") info; stop
  end if

! Orimo_ECS
  do irad = 1, dimbc
     v2jfc(irad) = (rrho2(irad, 1) + tmpr * bas_d2rpl1(irad, 0)) * bas_d2fac2(irad, 0) &
                 + (rrho2(irad, 2) + tmpi * bas_d2rpl1(irad, 0)) * bas_d2fac2(irad, 0) * iunit
  end do
  charge = (tmpr + tmpi * iunit) / bas_d2fac1(0)
  do irad = dimbc + 1, dim
     !! mean field in the scaled region
     v2jfc(irad) = bas_d2crpl1(irad, 0) * charge
  end do
! Orimo_ECS
  deallocate(rrho2)

end subroutine hprod_mkv2jfc_ecs
!######################################################################
subroutine hprod_v2xfc_opt1(orbg, v2xfc)

  use, intrinsic :: iso_c_binding
  use mod_sph, only : nlat, wlat
  use mod_rad, only : nradfc, nrad, xrad
  use mod_const, only : czero
  use mod_ormas, only : nfcore2, nfcore, nfun
  use mod_hprod, only : rho2, v2ang, v2sph

  implicit none
  complex(c_double_complex), intent(in) :: orbg(1:(nrad-1), 1:nlat, 1:nfun)
  complex(c_double_complex), intent(inout) :: v2xfc(1:(nrad-1))

  integer(c_int) :: llfun, ulfun
  integer(c_int), external :: util_omp_iproc
  integer(c_int) :: iproc, lll, ull, ifun, jfun, ilat, irad
  complex(c_double_complex) :: optfac
  complex(c_double_complex) :: eigx, eigy, eigxp, eigyp
  complex(c_double_complex), allocatable :: xwfn(:,:,:), ywfn(:,:,:)

  llfun = nfcore + 1
  ulfun = nfun

  allocate(xwfn(1:nradfc, 1:nlat, llfun:ulfun)) ! exact exchange times orbital
  allocate(ywfn(1:nradfc, 1:nlat, llfun:ulfun)) ! local exchange times orbital

  call hprod_mkrho2p('fc2dyn', orbg, v2ang)
  call bas_ang2sph2_fc2dyn(v2ang, rho2);
  call hprod_mkv2mf_poisson('fc2dyn', rho2, v2sph)
  call bas_sph2ang2_fc2dyn(v2sph, v2ang);
  call hprod_mkv2mf_herm('fc2dyn', v2sph, v2ang)

  eigx = czero
  eigy = czero
  !$omp parallel default(shared) private(iproc, eigxp, eigyp, lll, ull) reduction(+:eigx, eigy)
  !###########################
  iproc = util_omp_iproc()
  call util_omp_disp(1, nlat, lll, ull)

  eigxp = czero
  eigyp = czero
  xwfn(1:nradfc, lll:ull, llfun:ulfun) = czero
  ywfn(1:nradfc, lll:ull, llfun:ulfun) = czero
  do ifun = llfun, ulfun
     ! exact exchange
     do jfun = 1, nfcore2
        do ilat = lll, ull
           do irad = 1, nradfc
              xwfn(irad, ilat, ifun) = xwfn(irad, ilat, ifun) - v2ang(irad, ilat, jfun, ifun) * orbg(irad, ilat, jfun)
           end do
        end do
     end do

     ! local exchange
     do ilat = lll, ull
        do irad = 1, nradfc
           ywfn(irad, ilat, ifun) = ywfn(irad, ilat, ifun) + v2xfc(irad) * orbg(irad, ilat, ifun)
        end do
     end do

     ! contribution to orbital energy
     do ilat = lll, ull
        do irad = 1, nradfc
           eigxp = eigxp + conjg(orbg(irad, ilat, ifun)) * xwfn(irad, ilat, ifun) * wlat(ilat)
           eigyp = eigyp + conjg(orbg(irad, ilat, ifun)) * ywfn(irad, ilat, ifun) * wlat(ilat)
        end do
     end do
  end do

  eigx = eigx + eigxp
  eigy = eigy + eigyp
  !###########################
  !$omp end parallel

  optfac = eigx / eigy
  v2xfc(1:(nrad-1)) = v2xfc(1:(nrad-1)) * optfac

  !debug
  write(6, "('hprod_v2xfc_opt1: optfac = ', 6f20.10)") eigx, eigy, optfac
  do ifun = llfun, ulfun
     write(6, "('# ifun = ', i5)") ifun
     do irad = 1, nradfc
        write(6, "(f10.5)", advance = 'no') xrad(irad)
        do ilat = 1, 1
           write(6, "(6e20.10)", advance = 'no') &
                xwfn(irad, ilat, ifun), &
                optfac*ywfn(irad, ilat, ifun), &
                optfac*ywfn(irad, ilat, ifun)-xwfn(irad, ilat, ifun)
        end do
        write(6, *)
     end do
     write(6, *)
     write(6, *)
  end do
  !debug

  deallocate(ywfn)
  deallocate(xwfn)

end subroutine hprod_v2xfc_opt1
!######################################################################
subroutine hprod_v2xfc_opt2(cic, orbg, v2xfc)

  use, intrinsic :: iso_c_binding
  use mod_sph, only : nlat, wlat
  use mod_rad, only : nradfc, nrad, xrad
  use mod_const, only : czero, ctwo
  use mod_hprod, only : rho2, v2ang, v2sph, den1, den2
  use mod_ormas, only : nfcore2, nfcore, ncore, nact, nfun

  implicit none
  complex(c_double_complex), intent(in) :: cic(*)
  complex(c_double_complex), intent(in) :: orbg(1:(nrad-1), 1:nlat, 1:nfun)
  complex(c_double_complex), intent(inout) :: v2xfc(1:(nrad-1))

  integer(c_int) :: llfun, ulfun
  integer(c_int), external :: util_omp_iproc
  integer(c_int) :: iproc, lll, ull, ilat, irad
  integer(c_int) :: ifun, jfun, kfun, lfun
  integer(c_int) :: iact, jact, kact, lact
  complex(c_double_complex) :: optfac, klkl, klkx
  complex(c_double_complex), allocatable :: xwfn(:,:,:), ywfn(:,:,:)
  complex(c_double_complex), allocatable :: kintx(:,:), kintl(:,:)
  complex(c_double_complex), allocatable :: int1f(:,:), int2f(:,:,:,:)
  complex(c_double_complex), allocatable :: int1a(:,:), int2a(:,:,:,:)
  real(c_double), external :: hprod_ene_act

  llfun = nfcore + 1
  ulfun = nfun

  allocate(xwfn(1:nradfc, 1:nlat, llfun:ulfun)) ! exact exchange times orbital
  allocate(ywfn(1:nradfc, 1:nlat, llfun:ulfun)) ! local exchange times orbital
  allocate(kintx(1:nfun, 1:nfun))
  allocate(kintl(1:nfun, 1:nfun))
  allocate(int1f(1:nfun, 1:nfun))
  allocate(int1a(1:nact, 1:nact))
  allocate(int2f(1:nfun, 1:nfun, 1:nfun, 1:nfun))
  allocate(int2a(1:nact, 1:nact, 1:nact, 1:nact))

  call ormas_mkden1(cic, den1)
  call ormas_mkden2(cic, den1, den2)

  call hprod_mkrho2p('fc2dyn', orbg, v2ang)
  call bas_ang2sph2_fc2dyn(v2ang, rho2);
  call hprod_mkv2mf_poisson('fc2dyn', rho2, v2sph)
  call bas_sph2ang2_fc2dyn(v2sph, v2ang);
  call hprod_mkv2mf_herm('fc2dyn', v2sph, v2ang)

  lll = 1
  ull = nlat

  xwfn(1:nradfc, lll:ull, llfun:ulfun) = czero
  ywfn(1:nradfc, lll:ull, llfun:ulfun) = czero

  do ifun = llfun, ulfun
     ! exact exchange
     do kfun = 1, nfcore2
        do ilat = lll, ull
           do irad = 1, nradfc
              xwfn(irad, ilat, ifun) = xwfn(irad, ilat, ifun) - v2ang(irad, ilat, kfun, ifun) * orbg(irad, ilat, kfun)
           end do
        end do
     end do
     ! local exchange
     do ilat = lll, ull
        do irad = 1, nradfc
           ywfn(irad, ilat, ifun) = ywfn(irad, ilat, ifun) + v2xfc(irad) * orbg(irad, ilat, ifun)
        end do
     end do
  end do

  ! <i|K_exx|j>
  kintx(1:nfun, 1:nfun) = czero
  do ifun = llfun, ulfun
     do jfun = llfun, ulfun
        do ilat = lll, ull
           do irad = 1, nradfc
              kintx(ifun, jfun) = kintx(ifun, jfun) + conjg(orbg(irad, ilat, ifun)) * xwfn(irad, ilat, jfun) * wlat(ilat)
           end do
        end do
     end do
  end do
  ! <i|K_loc|j>
  kintl(1:nfun, 1:nfun) = czero
  do ifun = llfun, ulfun
     do jfun = llfun, ulfun
        do ilat = lll, ull
           do irad = 1, nradfc
              kintl(ifun, jfun) = kintl(ifun, jfun) + conjg(orbg(irad, ilat, ifun)) * ywfn(irad, ilat, jfun) * wlat(ilat)
           end do
        end do
     end do
  end do

  ! ##### <Psi|K_loc*K_loc|Psi> #####
  ! <i|K_loc*K_loc|j>
  int1f(1:nfun, 1:nfun) = czero
  do ifun = llfun, ulfun
     do jfun = llfun, ulfun
        do ilat = lll, ull
           do irad = 1, nradfc
              int1f(ifun, jfun) = int1f(ifun, jfun) + conjg(ywfn(irad, ilat, ifun)) * ywfn(irad, ilat, jfun) * wlat(ilat)
           end do
        end do
     end do
  end do
  ! <i|K_loc|j><k|K_loc|l>
  do ifun = llfun, ulfun
     do jfun = llfun, ulfun
        do kfun = llfun, ulfun
           do lfun = llfun, ulfun
              int2f(ifun, jfun, kfun, lfun) = kintl(ifun, jfun) * kintl(kfun, lfun) * ctwo
           end do
        end do
     end do
  end do

  klkl = czero
  ! dynamical-core contributions
  do ifun = nfcore + 1, ncore
     klkl = klkl + int1f(ifun, ifun) * ctwo
     do jfun = nfcore + 1, ncore
        klkl = klkl + int2f(ifun, ifun, jfun, jfun) * ctwo - int2f(ifun, jfun, jfun, ifun)
     end do
  end do
  ! active contributions
  int1a(1:nact, 1:nact) = czero
  do iact = 1, nact
     ifun = ncore + iact
     do jact = 1, nact
        jfun = ncore + jact
        int1a(iact, jact) = int1f(ifun, jfun)
        do kfun = nfcore + 1, ncore
           int1a(iact, jact) = int1a(ifun, jfun) + int2f(ifun, jfun, kfun, kfun) * ctwo - int2f(ifun, kfun, kfun, jfun)
        end do
        do kact = 1, nact
           kfun = ncore + kact
           do lact = 1, nact
              lfun = ncore + lact
              int2a(iact, jact, kact, lact) = int2f(ifun, jfun, kfun, lfun)
           end do
        end do
     end do
  end do
  klkl = klkl + hprod_ene_act(int1a, int2a, den1, den2)

  ! ##### <Psi|K_loc*K_exx|Psi> #####
  ! <i|K_loc*K_exx|j>
  int1f(1:nfun, 1:nfun) = czero
  do ifun = llfun, ulfun
     do jfun = llfun, ulfun
        do ilat = lll, ull
           do irad = 1, nradfc
              int1f(ifun, jfun) = int1f(ifun, jfun) + conjg(ywfn(irad, ilat, ifun)) * xwfn(irad, ilat, jfun) * wlat(ilat)
           end do
        end do
     end do
  end do
  ! <i|K_loc|j><k|K_exx|l>
  do ifun = llfun, ulfun
     do jfun = llfun, ulfun
        do kfun = llfun, ulfun
           do lfun = llfun, ulfun
              int2f(ifun, jfun, kfun, lfun) = kintl(ifun, jfun) * kintx(kfun, lfun) * ctwo
           end do
        end do
     end do
  end do

  klkx = czero
  ! dynamical-core contributions
  do ifun = nfcore + 1, ncore
     klkx = klkx + int1f(ifun, ifun) * ctwo
     do jfun = nfcore + 1, ncore
        klkx = klkx + int2f(ifun, ifun, jfun, jfun) * ctwo - int2f(ifun, jfun, jfun, ifun)
     end do
  end do
  ! active contributions
  int1a(1:nact, 1:nact) = czero
  do iact = 1, nact
     ifun = ncore + iact
     do jact = 1, nact
        jfun = ncore + jact
        int1a(iact, jact) = int1f(ifun, jfun)
        do kfun = nfcore + 1, ncore
           int1a(iact, jact) = int1a(ifun, jfun) + int2f(ifun, jfun, kfun, kfun) * ctwo - int2f(ifun, kfun, kfun, jfun)
        end do
        do kact = 1, nact
           kfun = ncore + kact
           do lact = 1, nact
              lfun = ncore + lact
              int2a(iact, jact, kact, lact) = int2f(ifun, jfun, kfun, lfun)
           end do
        end do
     end do
  end do
  klkx = klkx + hprod_ene_act(int1a, int2a, den1, den2)

  optfac = klkx / klkl
  v2xfc(1:(nrad-1)) = v2xfc(1:(nrad-1)) * optfac

  !debug
  write(6, "('hprod_v2xfc_opt2: optfac = ', 6f20.10)") klkx, klkl, optfac
  do ifun = llfun, ulfun
     write(6, "('# ifun = ', i5)") ifun
     do irad = 1, nradfc
        write(6, "(f10.5)", advance = 'no') xrad(irad)
        do ilat = 1, 1
           write(6, "(6e20.10)", advance = 'no') &
                xwfn(irad, ilat, ifun), &
                optfac*ywfn(irad, ilat, ifun), &
                optfac*ywfn(irad, ilat, ifun)-xwfn(irad, ilat, ifun)
        end do
        write(6, *)
     end do
     write(6, *)
     write(6, *)
  end do
  !debug

  deallocate(int2a)
  deallocate(int2f)
  deallocate(int1a)
  deallocate(int1f)
  deallocate(kintl)
  deallocate(kintx)
  deallocate(ywfn)
  deallocate(xwfn)

end subroutine hprod_v2xfc_opt2
!######################################################################
subroutine hprod_v2fc_enefc(enefc)

  use, intrinsic :: iso_c_binding
  use mod_control, only : dft_type, exact3j
  use mod_sph, only : lmax1, nlat
  use mod_rad, only : nrad, nradfc
  use mod_ormas, only : nfcore, nfun
  use mod_const, only : czero, two, ctrue
  use mod_hprod, only : orb, orbg, h0orb, gorb, gorbg, v2ang, v2sph
!debug
  use mod_rad, only : xrad
  use mod_bas, only : lval, bas_d2crpl1
!debug

  implicit none
  real(c_double), intent(out) :: enefc
! real(c_double), external :: hprod_ene_fcx
  integer(c_int) :: ifun, irad, l
  complex(c_double_complex) :: tmp

  h0orb(1:nradfc, 0:lmax1, 1:nfcore) = czero
  call hprod_tprod_fc(orb, h0orb)

  tmp = czero
  do ifun = 1, nfcore
     do l = 0, lmax1
        do irad = 1, nradfc
           tmp = tmp + conjg(orb(irad, l, ifun)) * h0orb(irad, l, ifun)
        end do
     end do
  end do
  enefc = dble(tmp) * two

  gorbg(1:nradfc, 1:nlat, 1:nfcore) = czero
  if (.not.exact3j) then
     !call hprod_mfprod_fcj(ctrue, gorbg, 1, nlat)
     call hprod_mfprod2_fcj(ctrue, orbg, gorbg)
     if (dft_type == 0) then
        !call hprod_mfprod_fcx2(ctrue, gorbg, 1, nlat)
        !call hprod_mfprod_fcx1(ctrue, gorbg, 1, nlat)
        call hprod_mfprod2_fcx2(ctrue, orbg, gorbg)
        call hprod_mfprod2_fcx1(ctrue, orbg, v2ang, gorbg)
     end if
     call bas_ang2sph1_fc(gorbg, gorb);
  else
     call hprod_mfprodx3j_fcj(ctrue, orb, gorb)
     if (dft_type == 0) then
        call hprod_mfprodx3j_fcx2(ctrue, orb, gorb)
        call hprod_mfprodx3j_fcx1(ctrue, orb, v2sph, gorb)
     end if
  end if

  tmp = czero
  do ifun = 1, nfcore
     do l = 0, lmax1
        do irad = 1, nradfc
           tmp = tmp + conjg(orb(irad, l, ifun)) * gorb(irad, l, ifun)
        end do
     end do
  end do
  enefc = enefc + dble(tmp)
  write(6, "('# hprod_v2fc_enefc: enefc = ', f20.10)") enefc

!debug
!  write(6, "('# hprod_v2fc_enefc: gorb')")
!  do irad = 1, nrad - 1
!     write(6,"(f10.5)",advance='no') xrad(irad)
!     do ifun = 1, nfcore
!        write(6,"(6f20.10)",advance='no') gorb(irad,0,ifun),v2sph(irad,0,ifun,ifun),bas_d2crpl1(irad,0)
!     end do
!     write(6,*)
!  end do
!  stop
!debug

end subroutine hprod_v2fc_enefc
!######################################################################
subroutine hprod_v2fc_xkli(v2xfc)

  use, intrinsic :: iso_c_binding
  use mod_sph, only : lmax1, nlat
  use mod_bas, only : mval, lval
  use mod_rad, only : nradfc, nrad, xrad, wrad
  use mod_ormas, only : nfcore2, nfcore, nfun
  use mod_hprod, only : orb, orbg, gorb, gorbg, rho2, v2sph, v2ang
  use mod_const, only : two, czero, runit, four, pi, CTRUE

  implicit none
  complex(c_double_complex), intent(out) :: v2xfc(1:(nrad-1))

  integer(c_int) :: llfun, ulfun
  integer(c_int) :: lll, ull, ifun, jfun, ilat, irad, nvar, ivar, jvar
  real(c_double), parameter :: small = 1.D-15
  real(c_double), parameter :: fourpi = four * pi
  complex(c_double_complex) :: errsla, esum1_sla, esum2_sla
  complex(c_double_complex) :: errkli, esum1_kli, esum2_kli
  complex(c_double_complex) :: tmp, lfac, slawfn, kliwfn
  integer(c_int), allocatable :: ind_map(:)
  complex(c_double_complex), allocatable :: dtot(:), dens(:,:)
  complex(c_double_complex), allocatable :: xsla(:), xkli(:)
  complex(c_double_complex), allocatable :: eigsla(:), eigexx(:), eigkli(:), amat(:,:)
  complex(c_double_complex), allocatable :: exxwfn(:,:)

  if (nfcore2 == 0) return
  llfun = nfcore + 1
  ulfun = nfun

  allocate(ind_map(1:nfun))
  ind_map(1:nfun) = 0

  nvar = 0
  do ifun = llfun, ulfun
     if (mval(ifun) == 0) then
        nvar = nvar + 1
        ind_map(nvar) = ifun
     end if
  end do

  allocate(dtot(1:(nrad-1)))
  allocate(dens(1:(nrad-1), 1:nvar))
  allocate(xsla(1:(nrad-1)))
  allocate(xkli(1:(nrad-1)))
  allocate(eigsla(1:nvar))
  allocate(eigexx(1:nvar))
  allocate(eigkli(1:nvar))
  allocate(amat(1:nvar, 1:nvar))
  allocate(exxwfn(1:(nrad-1), 1:nvar))

  dtot(1:(nrad-1)) = czero
  dens(1:(nrad-1), 1:nvar) = czero
  xsla(1:(nrad-1)) = czero
  xkli(1:(nrad-1)) = czero
  eigsla(1:nvar) = czero
  eigexx(1:nvar) = czero
  eigkli(1:nvar) = czero
  amat(1:nvar, 1:nvar) = czero
  exxwfn(1:(nrad-1), 1:nvar) = czero

  ! exact exchange
  call hprod_mkrho2p('fc2dyn', orbg, v2ang)
  call bas_ang2sph2_fc2dyn(v2ang, rho2);
  call hprod_mkv2mf_poisson('fc2dyn', rho2, v2sph)
  call bas_sph2ang2_fc2dyn(v2sph, v2ang);
  call hprod_mkv2mf_herm('fc2dyn', v2sph, v2ang)
  gorbg(1:(nrad-1), 1:nlat, 1:nfun) = czero
  !$omp parallel default(shared) private(ifun, lll, ull)
  !###########################
  call util_omp_disp(1, nlat, lll, ull)
  do ifun = llfun, ulfun
     do jfun = 1, nfcore2
        do ilat = lll, ull
           do irad = 1, nradfc
              gorbg(irad, ilat, ifun) = gorbg(irad, ilat, ifun) - v2ang(irad, ilat, jfun, ifun) * orbg(irad, ilat, jfun)
           end do
        end do
     end do
  end do
  !###########################
  !$omp end parallel
  call bas_ang2sph1_dyn(gorbg, gorb);
  do ivar = 1, nvar
     ifun = ind_map(ivar)
     exxwfn(1:nradfc, ivar) = gorb(1:nradfc, lval(ifun), ifun)
  end do

!DEBUG
!  write(6, "('hprod_v2fc_xkli: exxwfn')")
!  do irad = 1, nradfc
!     write(6, "(i5, f12.5)", advance = 'no') irad, xrad(irad)
!     do ivar = 1, nvar
!        ifun = ind_map(ivar)
!        write(6, "(f12.5)", advance = 'no') dble(exxwfn(irad, ivar))
!     end do
!     write(6, *)
!  end do
!DEBUG

  ! density
  do ivar = 1, nvar
     ifun = ind_map(ivar)
     lfac = (2 * lval(ifun) + 1) / fourpi
     do irad = 1, nradfc
        tmp = conjg(orb(irad, lval(ifun), ifun)) * orb(irad, lval(ifun), ifun)
        dens(irad, ivar) = tmp
        dtot(irad) = dtot(irad) + tmp * lfac
     end do
  end do

!DEBUG
!  write(6, "('hprod_v2fc_xkli: dens')")
!  do irad = 1, nradfc
!     write(6, "(i5, f12.5)", advance = 'no') irad, xrad(irad)
!     do ivar = 1, nvar
!        ifun = ind_map(ivar)
!        write(6, "(f12.5)", advance = 'no') dble(dens(irad, ivar))
!     end do
!     write(6, "(f12.5)") dble(dtot(irad))
!  end do
!DEBUG

  ! matrix elements of exact exchange
  ! slater's average exchange
  do ivar = 1, nvar
     ifun = ind_map(ivar)
     lfac = (2 * lval(ifun) + 1) / fourpi
     do irad = 1, nradfc
        tmp = conjg(orb(irad, lval(ifun), ifun)) * exxwfn(irad, ivar)
        eigexx(ivar) = eigexx(ivar) + tmp
        xsla(irad) = xsla(irad) + lfac * tmp
     end do
  end do
  do irad = 1, nradfc
     if (abs(dtot(irad)) > small) then
        xsla(irad) = xsla(irad) / dtot(irad)
     end if
  end do

  ! matrix elements of slater's average exchange
  do ivar = 1, nvar
     do irad = 1, nradfc
        eigsla(ivar) = eigsla(ivar) + xsla(irad) * dens(irad, ivar)
     end do
  end do

  ! cross correlation matrix
  do ivar = 1, nvar
     do jvar = 1, nvar
        jfun = ind_map(jvar)
        lfac = (2 * lval(jfun) + 1) / fourpi
        do irad = 1, nradfc
           if (abs(dtot(irad)) > small) then
              amat(ivar, jvar) = amat(ivar, jvar) - dens(irad, ivar) * dens(irad, jvar) / dtot(irad) * lfac
           end if
        end do
     end do
  end do

  ! right-hand side
  do ivar = 1, nvar
     eigkli(ivar) = eigsla(ivar)
     do jvar = 1, nvar
        eigkli(ivar) = eigkli(ivar) + amat(ivar, jvar) * eigexx(jvar)
     end do
  end do

  ! coefficient matrix
  do ivar = 1, nvar
     amat(ivar, ivar) = runit + amat(ivar, ivar)
  end do
  
!debug
  write(6, "('amat & bvec:')")
  do ivar = 1, nvar
     write(6, "(2i5)", advance = 'no') ivar, ind_map(ivar)
     do jvar = 1, nvar
        write(6, "(e20.10)", advance = 'no') dble(amat(ivar, jvar))
     end do
     write(6, "(' | ')", advance = 'no')
     write(6, "(e20.10)") dble(eigkli(ivar))
  end do
!debug

!DEBUG
!  call lapack_zheev(nvar, amat, umat)
!  write(6, "('eigensolutions of amat:')")
!  do ivar = 1, nvar
!     write(6, "(i5)", advance = 'no') ivar
!     write(6, "(e20.10, ' | ')", advance = 'no') dble(amat(ivar, ivar))
!     do jvar = 1, nvar
!        write(6, "(e20.10)", advance = 'no') dble(umat(jvar, ivar))
!     end do
!     write(6, *)
!  end do  
!  stop
!DEBUG

  ! solve linear equation
!  call futil_lineq_svd(nvar, 1.D-10, amat, eigkli)
!  call lapack_zhesv(nvar, 1, amat, eigkli)
  call lapack_zgesv(nvar, 1, amat, eigkli)

!debug
  write(6, "('orbital exchange energy:')")
  do ivar = 1, nvar
     write(6, "(2i5, 6e20.10)") ivar, ind_map(ivar), eigsla(ivar), eigkli(ivar), eigexx(ivar)
  end do
!debug

  ! kli exchange
  xkli(1:nradfc) = xsla(1:nradfc)
  if (.false.) then
     write(6, "('hprod_v2fc_xkli: KLI may have bug...')")
     do ivar = 1, nvar
        ifun = ind_map(ivar)
        lfac = (2 * lval(ifun) + 1) / fourpi
        do irad = 1, nradfc
           if (abs(dtot(irad)) > small) then
              xkli(irad) = xkli(irad) + dens(irad, ivar) * (eigkli(ivar) - eigexx(ivar)) / dtot(irad) * lfac
           end if
        end do
     end do
  end if

!debug
  write(6, "('orbital exchange energy:')")
  do ivar = 1, nvar
     tmp = czero
     do irad = 1, nradfc
        tmp = tmp + dens(irad, ivar) * xkli(irad)
     end do
     write(6, "(2i5, 6e20.10)") ivar, ind_map(ivar), eigsla(ivar), tmp, eigexx(ivar)
  end do

  write(6, "('# hprod_v2fc_xkli: local-fcx')")
  do irad = 1, nradfc
     write(6, "(3e20.10)") xrad(irad), dble(xsla(irad)), dble(xkli(irad))
  end do
  write(6, *)
  write(6, *)

  write(6, "('# hprod_v2fc_xkli: exchange times orbitals')")
  esum2_sla = czero
  esum2_kli = czero
  do ifun = llfun, ulfun
     esum1_sla = czero
     esum1_kli = czero
     write(6, "('# ifun = ', i5)") ifun
     do irad = 1, nradfc
        slawfn = xsla(irad) * orb(irad, lval(ifun), ifun)
        kliwfn = xkli(irad) * orb(irad, lval(ifun), ifun)
        write(6, "(5e20.10)") xrad(irad), dble(orb(irad, lval(ifun), ifun)), &
             dble(slawfn), dble(kliwfn), dble(gorb(irad, lval(ifun), ifun))
        errsla = slawfn - gorb(irad, lval(ifun), ifun)
        errkli = kliwfn - gorb(irad, lval(ifun), ifun)
        esum1_sla = esum1_sla + conjg(errsla) * errsla
        esum1_kli = esum1_kli + conjg(errkli) * errkli
     end do
     esum2_sla = esum2_sla + esum1_sla
     esum2_kli = esum2_kli + esum1_kli
     write(6, "('# esum1(',i5,') =', 2f20.10)") ifun, dble(esum1_sla), dble(esum1_kli)
     write(6, *)
     write(6, *)
  end do
  write(6, "('# esum2 =', 2f20.10)") dble(esum2_sla), dble(esum2_kli)
!debug

  v2xfc(1:(nrad-1)) = czero
  v2xfc(1:nradfc) = xkli(1:nradfc)

  deallocate(ind_map)
  deallocate(dtot)
  deallocate(dens)
  deallocate(xsla)
  deallocate(xkli)
  deallocate(eigsla)
  deallocate(eigexx)
  deallocate(eigkli)
  deallocate(amat)
  deallocate(exxwfn)
  return

end subroutine hprod_v2fc_xkli
!#######################################################################
