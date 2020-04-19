!######################################################################
subroutine hprod_mkrho2()

  use, intrinsic :: iso_c_binding
  use mod_hprod, only : orb, orbg, rho2, v2ang

  implicit none
  call hprod_mkrho2p('tot', orbg, v2ang)
  call bas_ang2sph2(v2ang, rho2);

end subroutine hprod_mkrho2
!######################################################################
subroutine hprod_mkrho2_dyn()

  use, intrinsic :: iso_c_binding
  use mod_hprod, only : orb, orbg, rho2, v2ang
  use mod_hprod, only : orbe, orbo, v2ange, v2ango, rho23j
!debug
!  use mod_ormas, only : nfun
!  use mod_rad, only : nrad, wrad
!  use mod_sph, only : lmax2, nlat, wlat
!  use mod_bas, only : mval
!debug

  implicit none
!debug
!  integer(c_long) :: ifun,jfun,l,ilat,irad
!  complex(c_double_complex), allocatable :: test_ovlp(:,:)
!debug

  call hprod_mkrho2p('dyn', orbg, v2ang)
  call bas_ang2sph2_dyn(v2ang, rho2);
  call hprod_mkrho2p('fc1dyn', orbg, v2ang)
  call bas_ang2sph2_fc1dyn(v2ang, rho2);

!debug
!  write(6, "('hprod_mkrho2_dyn: rho2')")
!  call hprod_printrho(rho2)
!  stop
!
!  allocate(test_ovlp(1:nfun,1:nfun))
!  call hprod_mkovlp(1D10, orb, orb, test_ovlp)
!  write(6, "('ovlp from orbsph:')")
!  do ifun = 1, nfun
!  do jfun = 1, nfun
!     write(6,"(2i5,2f20.10)") ifun,jfun,test_ovlp(ifun,jfun)
!  end do
!  end do
!
!  test_ovlp = 0d0
!  do ifun = 1, nfun
!  do jfun = 1, nfun
!     do ilat = 1, nlat
!     do irad = 1, nrad-1
!        test_ovlp(ifun,jfun) = test_ovlp(ifun,jfun) + conjg(orbg(irad,ilat,ifun))*orbg(irad,ilat,jfun)*wlat(ilat)
!     end do
!     end do
!  end do
!  end do
!  write(6, "('ovlp from orbang:')")
!  do ifun = 1, nfun
!  do jfun = 1, nfun
!     write(6,"(2i5,2f20.10)") ifun,jfun,test_ovlp(ifun,jfun)
!  end do
!  end do
!
!  test_ovlp = 0d0
!  do ifun = 1, nfun
!  do jfun = 1, nfun
!     do ilat = 1, nlat
!     do irad = 1, nrad-1
!        test_ovlp(ifun,jfun) = test_ovlp(ifun,jfun) + v2ang(irad,ilat,ifun,jfun)*wlat(ilat)
!     end do
!     end do
!  end do
!  end do
!  write(6, "('ovlp from rhoang:')")
!  do ifun = 1, nfun
!  do jfun = 1, nfun
!     write(6,"(2i5,2f20.10)") ifun,jfun,test_ovlp(ifun,jfun)
!  end do
!  end do
!
!  test_ovlp = 0d0
!  do ifun = 1, nfun
!  do jfun = 1, nfun
!     do l = 0, lmax2
!     do irad = 1, nrad-1
!        test_ovlp(ifun,jfun) = test_ovlp(ifun,jfun) + rho2(irad,l,ifun,jfun)
!     end do
!     end do
!  end do
!  end do
!  write(6, "('ovlp from rhosph:')")
!  do ifun = 1, nfun
!  do jfun = 1, nfun
!     write(6,"(2i5,2f20.10)") ifun,jfun,test_ovlp(ifun,jfun)
!  end do
!  end do
!
!  deallocate(test_ovlp)
!debug

! ##### 3j selection rule #####
!old  if (exact3j) then
!old     call hprod_mkrho2p3j('dyn', orbe, orbo, v2ange, v2ango)
!old     call bas_ang2sph2_dyn3j(v2ange, v2ango, rho23j);
!old     call hprod_mkrho2p3j('fc1dyn', orbe, orbo, v2ange, v2ango)
!old     call bas_ang2sph2_fc1dyn3j(v2ange, v2ango, rho23j);
!old  end if
! ##### 3j selection rule #####

end subroutine hprod_mkrho2_dyn
!######################################################################
subroutine hprod_mkrho2p(v2_type, orbg, rho2g)

  use, intrinsic :: iso_c_binding
  use mod_sph, only : nlat
  use mod_const, only : ctwo
  use mod_ormas, only : nfun, nfcore, nfcore1, nfcore2
  use mod_rad, only : xrad, wrad, nrad, nradfc

  implicit none
  character(len=*), intent(in) :: v2_type
  complex(c_double_complex), intent(in) :: orbg(1:(nrad-1), 1:nlat, 1:nfun)
  complex(c_double_complex), intent(out) :: rho2g(1:(nrad-1), 1:nlat, 1:nfun, 1:nfun)
  integer(c_long) :: dim, ifun_ll, ifun_ul, jfun_ll, jfun_ul, llr, ulr, ifun, jfun, ilat, irad

  if (trim(v2_type) == 'tot') then
     dim = nrad - 1
     ifun_ll = 1
     ifun_ul = nfun
     jfun_ll = 1
     jfun_ul = nfun
  else if (trim(v2_type) == 'fcfc') then
     dim = nradfc
!    dim = nrad - 1
     ifun_ll = 1
     ifun_ul = nfcore
     jfun_ll = 1
     jfun_ul = nfcore
  else if (trim(v2_type) == 'fc1fc1') then
     dim = nradfc
!    dim = nrad - 1
     ifun_ll = nfcore2 + 1
     ifun_ul = nfcore
     jfun_ll = nfcore2 + 1
     jfun_ul = nfcore
  else if (trim(v2_type) == 'fc2fc2') then
     dim = nradfc
!    dim = nrad - 1
     ifun_ll = 1
     ifun_ul = nfcore2
     jfun_ll = 1
     jfun_ul = nfcore2
  else if (trim(v2_type) == 'fc1dyn') then
     dim = nradfc
     ifun_ll = nfcore2 + 1
     ifun_ul = nfcore
     jfun_ll = nfcore + 1
     jfun_ul = nfun
  else if (trim(v2_type) == 'fc2dyn') then
     dim = nradfc
     ifun_ll = 1
     ifun_ul = nfcore2
     jfun_ll = nfcore + 1
     jfun_ul = nfun
  else if (trim(v2_type) == 'dyn') then
     dim = nrad - 1
     ifun_ll = nfcore + 1
     ifun_ul = nfun
     jfun_ll = nfcore + 1
     jfun_ul = nfun
  else
     stop "bad v2_type in hprod_mkrho2p"
  end if

  !$omp parallel default(shared) private(llr, ulr)
  !###########################
  call util_omp_disp(1, dim, llr, ulr)
  do ifun = ifun_ll, ifun_ul
     do jfun = max(jfun_ll, ifun), max(jfun_ul, ifun)
!     do jfun = jfun_ll, jfun_ul
        do ilat = 1, nlat
           do irad = llr, ulr
             rho2g(irad, ilat, jfun, ifun) = conjg(orbg(irad, ilat, jfun)) &
                                                 * orbg(irad, ilat, ifun)
!            rho2g(irad, ilat, jfun, ifun) = conjg(orbg(irad, ilat, jfun)) &
!                                                * orbg(irad, ilat, ifun) * wrad(irad)
           end do
        end do
     end do
  end do
  !###########################
  !$omp end parallel

end subroutine hprod_mkrho2p
!######################################################################
subroutine hprod_mkrho2p3j(v2_type, orbe, orbo, rho2e, rho2o)

  use, intrinsic :: iso_c_binding
  use mod_sph, only : nlat
  use mod_const, only : ctwo
  use mod_ormas, only : nfun, nfcore, nfcore1, nfcore2
  use mod_rad, only : xrad, wrad, nrad, nradfc

  implicit none
  character(len=*), intent(in) :: v2_type
  complex(c_double_complex), intent(in) :: orbe(1:(nrad-1), 1:nlat, 1:nfun)
  complex(c_double_complex), intent(in) :: orbo(1:(nrad-1), 1:nlat, 1:nfun)
  complex(c_double_complex), intent(out) :: rho2e(1:(nrad-1), 1:nlat, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(out) :: rho2o(1:(nrad-1), 1:nlat, 1:nfun, 1:nfun)
  integer(c_long) :: dim, ifun_ll, ifun_ul, jfun_ll, jfun_ul, llr, ulr, ifun, jfun, ilat, irad

  if (trim(v2_type) == 'tot') then
     dim = nrad - 1
     ifun_ll = 1
     ifun_ul = nfun
     jfun_ll = 1
     jfun_ul = nfun
  else if (trim(v2_type) == 'fcfc') then
     dim = nradfc
!    dim = nrad - 1
     ifun_ll = 1
     ifun_ul = nfcore
     jfun_ll = 1
     jfun_ul = nfcore
  else if (trim(v2_type) == 'fc1fc1') then
     dim = nradfc
!    dim = nrad - 1
     ifun_ll = nfcore2 + 1
     ifun_ul = nfcore
     jfun_ll = nfcore2 + 1
     jfun_ul = nfcore
  else if (trim(v2_type) == 'fc2fc2') then
     dim = nradfc
!    dim = nrad - 1
     ifun_ll = 1
     ifun_ul = nfcore2
     jfun_ll = 1
     jfun_ul = nfcore2
  else if (trim(v2_type) == 'fc1dyn') then
     dim = nradfc
     ifun_ll = nfcore2 + 1
     ifun_ul = nfcore
     jfun_ll = nfcore + 1
     jfun_ul = nfun
  else if (trim(v2_type) == 'fc2dyn') then
     dim = nradfc
     ifun_ll = 1
     ifun_ul = nfcore2
     jfun_ll = nfcore + 1
     jfun_ul = nfun
  else if (trim(v2_type) == 'dyn') then
     dim = nrad - 1
     ifun_ll = nfcore + 1
     ifun_ul = nfun
     jfun_ll = nfcore + 1
     jfun_ul = nfun
  else
     stop "bad v2_type in hprod_mkrho2p3j"
  end if

  !$omp parallel default(shared) private(llr, ulr)
  !###########################
  call util_omp_disp(1, dim, llr, ulr)
  do ifun = ifun_ll, ifun_ul
     do jfun = max(jfun_ll, ifun), max(jfun_ul, ifun)
        do ilat = 1, nlat
           do irad = llr, ulr
             rho2e(irad, ilat, jfun, ifun) = conjg(orbe(irad, ilat, jfun)) &
                                                 * orbe(irad, ilat, ifun) &
                                           + conjg(orbo(irad, ilat, jfun)) &
                                                 * orbo(irad, ilat, ifun)
             rho2o(irad, ilat, jfun, ifun) = conjg(orbe(irad, ilat, jfun)) &
                                                 * orbo(irad, ilat, ifun) &
                                           + conjg(orbo(irad, ilat, jfun)) &
                                                 * orbe(irad, ilat, ifun)
           end do
        end do
     end do
  end do
  !###########################
  !$omp end parallel

end subroutine hprod_mkrho2p3j
!######################################################################
