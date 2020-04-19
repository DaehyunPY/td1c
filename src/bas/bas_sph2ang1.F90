!######################################################################
subroutine bas_sph2ang1(orb1, orb2)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad
  use mod_sph, only : lmax1, nlat, legf1
  use mod_bas, only : mval
  use mod_ormas, only : nfun
  use mod_const, only : czero

  implicit none
  complex(c_double_complex), intent(in) :: orb1(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(out) :: orb2(1:(nrad-1), 1:nlat, 1:nfun)

  integer(c_int) :: ifun, irad, ilat, l, mi, llr, ulr

  !$omp parallel default(shared) private(mi, llr, ulr)
  call util_omp_disp(1, nrad - 1, llr, ulr)

  do ifun = 1, nfun
     orb2(llr:ulr, 1:nlat, ifun) = czero
!    mi = abs(mval(ifun))
     mi = mval(ifun)
     do ilat = 1, nlat
        do l = abs(mi), lmax1
           do irad = llr, ulr
              orb2(irad, ilat, ifun) = orb2(irad, ilat, ifun) &
                                     + orb1(irad, l,    ifun) * legf1(l, ilat, mi)
           end do
        end do
     end do
!     ! Condon-Shortley factor
!     if (mval(ifun) > 0 .and. mod(mval(ifun),2) == 1) then
!        do ilat = 1, nlat
!           do irad = llr, ulr
!              orb2(irad, ilat, ifun) = - orb2(irad, ilat, ifun)
!           end do
!        end do
!     end if
  end do
  !$omp end parallel

end subroutine bas_sph2ang1
!######################################################################
subroutine bas_sph2ang1_radfc(orb1, orb2)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad, nradfc
  use mod_sph, only : lmax1, nlat, legf1
  use mod_bas, only : mval
  use mod_ormas, only : nfun
  use mod_const, only : czero

  implicit none
  complex(c_double_complex), intent(in) :: orb1(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(out) :: orb2(1:(nrad-1), 1:nlat, 1:nfun)

  integer(c_int) :: ifun, irad, ilat, l, mi, llr, ulr

  !$omp parallel default(shared) private(mi, llr, ulr)
  call util_omp_disp(1, nradfc, llr, ulr)

  do ifun = 1, nfun
     orb2(llr:ulr, 1:nlat, ifun) = czero
!    mi = abs(mval(ifun))
     mi = mval(ifun)
     do ilat = 1, nlat
        do l = abs(mi), lmax1
           do irad = llr, ulr
              orb2(irad, ilat, ifun) = orb2(irad, ilat, ifun) &
                                     + orb1(irad, l,    ifun) * legf1(l, ilat, mi)
           end do
        end do
     end do
!     ! Condon-Shortley factor
!     if (mval(ifun) > 0 .and. mod(mval(ifun),2) == 1) then
!        do ilat = 1, nlat
!           do irad = llr, ulr
!              orb2(irad, ilat, ifun) = - orb2(irad, ilat, ifun)
!           end do
!        end do
!     end if
  end do
  !$omp end parallel

end subroutine bas_sph2ang1_radfc
!######################################################################
subroutine bas_sph2ang1_dyn(orb1, orb2)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad
  use mod_sph, only : lmax1, nlat, legf1
  use mod_bas, only : mval
  use mod_ormas, only : nfcore, nfun
  use mod_const, only : czero

  implicit none
  complex(c_double_complex), intent(in) :: orb1(1:(nrad-1), 0:lmax1, 1:*)
  complex(c_double_complex), intent(out) :: orb2(1:(nrad-1), 1:nlat, 1:*)

  integer(c_int) :: ifun, irad, ilat, l, mi, llr, ulr
  integer(c_int) :: il, nil, mapil(1:2, 1:nfun*nlat)

  nil = 0
  do ifun = nfcore + 1, nfun
     do ilat = 1, nlat
        nil = nil + 1
        mapil(1, nil) = ifun
        mapil(2, nil) = ilat
     end do
  end do

  !$omp parallel default(shared) private(mi,ifun,ilat)
  !$omp do
  do il = 1, nil
     ifun = mapil(1, il)
     ilat = mapil(2, il)
     mi = mval(ifun)
     orb2(1:(nrad-1), ilat, ifun) = czero
     do l = abs(mi), lmax1
        do irad = 1, nrad - 1
           orb2(irad, ilat, ifun) = orb2(irad, ilat, ifun) &
                                  + orb1(irad, l,    ifun) * legf1(l, ilat, mi)
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel

!def  !$omp parallel default(shared) private(mi, llr, ulr)
!def  call util_omp_disp(1, nrad - 1, llr, ulr)
!def
!def  do ifun = nfcore + 1, nfun
!def     orb2(llr:ulr, 1:nlat, ifun) = czero
!def!    mi = abs(mval(ifun))
!def     mi = mval(ifun)
!def     do ilat = 1, nlat
!def        do l = abs(mi), lmax1
!def           do irad = llr, ulr
!def              orb2(irad, ilat, ifun) = orb2(irad, ilat, ifun) &
!def                                     + orb1(irad, l,    ifun) * legf1(l, ilat, mi)
!def           end do
!def        end do
!def     end do
!def!     ! Condon-Shortley factor
!def!     if (mval(ifun) > 0 .and. mod(mval(ifun),2) == 1) then
!def!        do ilat = 1, nlat
!def!           do irad = llr, ulr
!def!              orb2(irad, ilat, ifun) = - orb2(irad, ilat, ifun)
!def!           end do
!def!        end do
!def!     end if
!def  end do
!def  !$omp end parallel

end subroutine bas_sph2ang1_dyn
!######################################################################
subroutine bas_sph2ang1_dyn3j(orb1, orb2e, orb2o)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad
  use mod_sph, only : lmax1, nlat, legf1
  use mod_bas, only : mval
  use mod_ormas, only : nfcore, nfun
  use mod_const, only : czero

  implicit none
  complex(c_double_complex), intent(in) :: orb1(1:(nrad-1), 0:lmax1, 1:*)
  complex(c_double_complex), intent(out) :: orb2e(1:(nrad-1), 1:nlat, 1:*)
  complex(c_double_complex), intent(out) :: orb2o(1:(nrad-1), 1:nlat, 1:*)

  integer(c_int) :: ifun, irad, ilat, l, mi, llr, ulr
  integer(c_int) :: il, nil, mapil(1:2, 1:nfun*nlat)

  nil = 0
  do ifun = nfcore + 1, nfun
     do ilat = 1, nlat
        nil = nil + 1
        mapil(1, nil) = ifun
        mapil(2, nil) = ilat
     end do
  end do

  !$omp parallel default(shared) private(mi,ifun,ilat)
  !$omp do
  do il = 1, nil
     ifun = mapil(1, il)
     ilat = mapil(2, il)
     mi = mval(ifun)
     orb2e(1:(nrad-1), ilat, ifun) = czero
     orb2o(1:(nrad-1), ilat, ifun) = czero
     do l = abs(mi), lmax1
        if (mod(l,2) == 0) then
           do irad = 1, nrad - 1
              orb2e(irad, ilat, ifun) = orb2e(irad, ilat, ifun) &
                                      + orb1 (irad, l,    ifun) * legf1(l, ilat, mi)
           end do
        else
           do irad = 1, nrad - 1
              orb2o(irad, ilat, ifun) = orb2o(irad, ilat, ifun) &
                                      + orb1 (irad, l,    ifun) * legf1(l, ilat, mi)
           end do
        end if
     end do
  end do
  !$omp end do
  !$omp end parallel
end subroutine bas_sph2ang1_dyn3j
!######################################################################
subroutine bas_sph2ang1_dyn_radfc(orb1, orb2)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad, nradfc
  use mod_sph, only : lmax1, nlat, legf1
  use mod_bas, only : mval
  use mod_ormas, only : nfcore, nfun
  use mod_const, only : czero

  implicit none
  complex(c_double_complex), intent(in) :: orb1(1:(nrad-1), 0:lmax1, 1:*)
  complex(c_double_complex), intent(out) :: orb2(1:(nrad-1), 1:nlat, 1:*)

  integer(c_int) :: ifun, irad, ilat, l, mi, llr, ulr

  !$omp parallel default(shared) private(mi, llr, ulr)
  call util_omp_disp(1, nradfc, llr, ulr)

  do ifun = nfcore + 1, nfun
     orb2(llr:ulr, 1:nlat, ifun) = czero
!    mi = abs(mval(ifun))
     mi = mval(ifun)
     do ilat = 1, nlat
        do l = abs(mi), lmax1
           do irad = llr, ulr
              orb2(irad, ilat, ifun) = orb2(irad, ilat, ifun) &
                                     + orb1(irad, l,    ifun) * legf1(l, ilat, mi)
           end do
        end do
     end do
!     ! Condon-Shortley factor
!     if (mval(ifun) > 0 .and. mod(mval(ifun),2) == 1) then
!        do ilat = 1, nlat
!           do irad = llr, ulr
!              orb2(irad, ilat, ifun) = - orb2(irad, ilat, ifun)
!           end do
!        end do
!     end if
  end do
  !$omp end parallel

end subroutine bas_sph2ang1_dyn_radfc
!######################################################################
subroutine bas_sph2ang1_fc(orb1, orb2)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad, nradfc
  use mod_sph, only : lmax1, nlat, legf1
  use mod_bas, only : mval
  use mod_ormas, only : nfcore
  use mod_const, only : czero

  implicit none
  complex(c_double_complex), intent(in) :: orb1(1:(nrad-1), 0:lmax1, 1:*)
  complex(c_double_complex), intent(out) :: orb2(1:(nrad-1), 1:nlat, 1:*)

  integer(c_int) :: ifun, irad, ilat, l, mi, llr, ulr

  !$omp parallel default(shared) private(mi, llr, ulr)
  call util_omp_disp(1, nradfc, llr, ulr)

  do ifun = 1, nfcore
     orb2(llr:ulr, 1:nlat, ifun) = czero
!    mi = abs(mval(ifun))
     mi = mval(ifun)
     do ilat = 1, nlat
        do l = abs(mi), lmax1
           do irad = llr, ulr
              orb2(irad, ilat, ifun) = orb2(irad, ilat, ifun) &
                                     + orb1(irad, l,    ifun) * legf1(l, ilat, mi)
           end do
        end do
     end do
!     ! Condon-Shortley factor
!     if (mval(ifun) > 0 .and. mod(mval(ifun),2) == 1) then
!        do ilat = 1, nlat
!           do irad = llr, ulr
!              orb2(irad, ilat, ifun) = - orb2(irad, ilat, ifun)
!           end do
!        end do
!     end if
  end do
  !$omp end parallel

end subroutine bas_sph2ang1_fc
!######################################################################
subroutine bas_sph2ang1_fc3j(orb1, orb2e, orb2o)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad, nradfc
  use mod_sph, only : lmax1, nlat, legf1
  use mod_bas, only : mval
  use mod_ormas, only : nfcore
  use mod_const, only : czero

  implicit none
  complex(c_double_complex), intent(in) :: orb1(1:(nrad-1), 0:lmax1, 1:*)
  complex(c_double_complex), intent(out) :: orb2e(1:(nrad-1), 1:nlat, 1:*)
  complex(c_double_complex), intent(out) :: orb2o(1:(nrad-1), 1:nlat, 1:*)

  integer(c_int) :: ifun, irad, ilat, l, llr, ulr

  !$omp parallel default(shared) private(llr, ulr)
  call util_omp_disp(1, nradfc, llr, ulr)

  do ifun = 1, nfcore
     orb2e(llr:ulr, 1:nlat, ifun) = czero
     orb2o(llr:ulr, 1:nlat, ifun) = czero
     do ilat = 1, nlat
        do l = abs(mval(ifun)), lmax1
           if (mod(l,2) == 0) then
              do irad = llr, ulr
                 orb2e(irad, ilat, ifun) = orb2e(irad, ilat, ifun) &
                                         + orb1 (irad, l,    ifun) * legf1(l, ilat, mval(ifun))
              end do
           else
              do irad = llr, ulr
                 orb2o(irad, ilat, ifun) = orb2o(irad, ilat, ifun) &
                                         + orb1 (irad, l,    ifun) * legf1(l, ilat, mval(ifun))
              end do
           end if
        end do
     end do
  end do
  !$omp end parallel

end subroutine bas_sph2ang1_fc3j
!######################################################################
