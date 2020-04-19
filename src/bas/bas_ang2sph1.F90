!######################################################################
subroutine bas_ang2sph1(orb1, orb2)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad
  use mod_sph, only : lmax1, nlat, legb1
  use mod_bas, only : mval
  use mod_ormas, only : nfun
  use mod_const, only : czero

  implicit none
  complex(c_double_complex), intent(in) :: orb1(1:(nrad-1), 1:nlat, 1:nfun)
  complex(c_double_complex), intent(out) :: orb2(1:(nrad-1), 0:lmax1, 1:nfun)

  integer(c_long) :: ifun, irad, ilat, l, mi, llr, ulr

  !$omp parallel default(shared) private(mi, llr, ulr)
  call util_omp_disp(1, nrad - 1, llr, ulr)

  do ifun = 1, nfun
     orb2(llr:ulr, 0:lmax1, ifun) = czero
!    mi = abs(mval(ifun))
     mi = mval(ifun)
     do l = abs(mi), lmax1
        do ilat = 1, nlat
           do irad = llr, ulr
              orb2(irad, l, ifun) = orb2(irad, l,    ifun) &
                                  + orb1(irad, ilat, ifun) * legb1(ilat, l, mi)
           end do
        end do
     end do
!     ! Condon-Shortley factor
!     if (mval(ifun) > 0 .and. mod(mval(ifun),2) == 1) then
!        do l = abs(mi), lmax1
!           do irad = llr, ulr
!              orb2(irad, l, ifun) = - orb2(irad, l, ifun)
!           end do
!        end do
!     end if
  end do
  !$omp end parallel

end subroutine bas_ang2sph1
!######################################################################
subroutine bas_ang2sph1_radfc(orb1, orb2)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad, nradfc
  use mod_sph, only : lmax1, nlat, legb1
  use mod_bas, only : mval
  use mod_ormas, only : nfun
  use mod_const, only : czero

  implicit none
  complex(c_double_complex), intent(in) :: orb1(1:(nrad-1), 1:nlat, 1:nfun)
  complex(c_double_complex), intent(out) :: orb2(1:(nrad-1), 0:lmax1, 1:nfun)

  integer(c_long) :: ifun, irad, ilat, l, mi, llr, ulr

  !$omp parallel default(shared) private(mi, llr, ulr)
  call util_omp_disp(1, nradfc, llr, ulr)

  do ifun = 1, nfun
     orb2(llr:ulr, 0:lmax1, ifun) = czero
!    mi = abs(mval(ifun))
     mi = mval(ifun)
     do l = abs(mi), lmax1
        do ilat = 1, nlat
           do irad = llr, ulr
              orb2(irad, l, ifun) = orb2(irad, l,    ifun) &
                                  + orb1(irad, ilat, ifun) * legb1(ilat, l, mi)
           end do
        end do
     end do
!     ! Condon-Shortley factor
!     if (mval(ifun) > 0 .and. mod(mval(ifun),2) == 1) then
!        do l = abs(mi), lmax1
!           do irad = llr, ulr
!              orb2(irad, l, ifun) = - orb2(irad, l, ifun)
!           end do
!        end do
!     end if
  end do
  !$omp end parallel

end subroutine bas_ang2sph1_radfc
!######################################################################
subroutine bas_ang2sph1_dyn(orb1, orb2)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad
  use mod_sph, only : lmax1, nlat, legb1
  use mod_bas, only : mval
  use mod_ormas, only : nfcore, nfun
  use mod_const, only : czero

  implicit none
  complex(c_double_complex), intent(in) :: orb1(1:(nrad-1), 1:nlat, 1:*)
  complex(c_double_complex), intent(out) :: orb2(1:(nrad-1), 0:lmax1, 1:*)

  integer(c_long) :: ifun, irad, ilat, l, mi, llr, ulr
  integer(c_long) :: il, nil, mapil(1:2, 1:nfun*(lmax1+1))

  nil = 0
  do ifun = nfcore + 1, nfun
     mi = mval(ifun)
     do l = abs(mi), lmax1
        nil = nil + 1
        mapil(1, nil) = ifun
        mapil(2, nil) = l
     end do
  end do

  !$omp parallel default(shared) private(mi,ifun,l)
  !$omp do
  do il = 1, nil
     ifun = mapil(1, il)
     l    = mapil(2, il)
     mi = mval(ifun)
     orb2(1:(nrad-1), l, ifun) = czero
     do ilat = 1, nlat
        do irad = 1, nrad - 1
           orb2(irad, l, ifun) = orb2(irad, l,    ifun) &
                               + orb1(irad, ilat, ifun) * legb1(ilat, l, mi)
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel

!def  !$omp parallel default(shared) private(mi, llr, ulr)
!def  call util_omp_disp(1, nrad - 1, llr, ulr)
!def
!def  do ifun = nfcore + 1, nfun
!def     orb2(llr:ulr, 0:lmax1, ifun) = czero
!def!    mi = abs(mval(ifun))
!def     mi = mval(ifun)
!def     do l = abs(mi), lmax1
!def        do ilat = 1, nlat
!def           do irad = llr, ulr
!def              orb2(irad, l, ifun) = orb2(irad, l,    ifun) &
!def                                  + orb1(irad, ilat, ifun) * legb1(ilat, l, mi)
!def           end do
!def        end do
!def     end do
!def!     ! Condon-Shortley factor
!def!     if (mval(ifun) > 0 .and. mod(mval(ifun),2) == 1) then
!def!        do l = abs(mi), lmax1
!def!           do irad = llr, ulr
!def!              orb2(irad, l, ifun) = - orb2(irad, l, ifun)
!def!           end do
!def!        end do
!def!     end if
!def  end do
!def  !$omp end parallel

end subroutine bas_ang2sph1_dyn
!######################################################################
subroutine bas_ang2sph1_dyn3j(orb1e, orb1o, orb2)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad
  use mod_sph, only : lmax1, nlat, legb1
  use mod_bas, only : mval
  use mod_ormas, only : nfcore, nfun
  use mod_const, only : czero

  implicit none
  complex(c_double_complex), intent(in) :: orb1e(1:(nrad-1), 1:nlat, 1:*)
  complex(c_double_complex), intent(in) :: orb1o(1:(nrad-1), 1:nlat, 1:*)
  complex(c_double_complex), intent(out) :: orb2(1:(nrad-1), 0:lmax1, 1:*)

  integer(c_long) :: ifun, irad, ilat, l, mi, llr, ulr
  integer(c_long) :: il, nil, mapil(1:2, 1:nfun*(lmax1+1))

  nil = 0
  do ifun = nfcore + 1, nfun
     mi = mval(ifun)
     do l = abs(mi), lmax1
        nil = nil + 1
        mapil(1, nil) = ifun
        mapil(2, nil) = l
     end do
  end do

  !$omp parallel default(shared) private(mi,ifun,l)
  !$omp do
  do il = 1, nil
     ifun = mapil(1, il)
     l    = mapil(2, il)
     mi = mval(ifun)
     orb2(1:(nrad-1), l, ifun) = czero
     if (mod(l,2) == 0) then
        do ilat = 1, nlat
           do irad = 1, nrad - 1
              orb2(irad, l, ifun) = orb2 (irad, l,    ifun) &
                                  + orb1e(irad, ilat, ifun) * legb1(ilat, l, mi)
           end do
        end do
     else
        do ilat = 1, nlat
           do irad = 1, nrad - 1
              orb2(irad, l, ifun) = orb2 (irad, l,    ifun) &
                                  + orb1o(irad, ilat, ifun) * legb1(ilat, l, mi)
           end do
        end do
     end if
  end do
  !$omp end do
  !$omp end parallel

end subroutine bas_ang2sph1_dyn3j
!######################################################################
subroutine bas_ang2sph1_dyn_radfc(orb1, orb2)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad, nradfc
  use mod_sph, only : lmax1, nlat, legb1
  use mod_bas, only : mval
  use mod_ormas, only : nfcore, nfun
  use mod_const, only : czero

  implicit none
  complex(c_double_complex), intent(in) :: orb1(1:(nrad-1), 1:nlat, 1:*)
  complex(c_double_complex), intent(out) :: orb2(1:(nrad-1), 0:lmax1, 1:*)

  integer(c_long) :: ifun, irad, ilat, l, mi, llr, ulr

  !$omp parallel default(shared) private(mi, llr, ulr)
  call util_omp_disp(1, nradfc, llr, ulr)

  do ifun = nfcore + 1, nfun
     orb2(llr:ulr, 0:lmax1, ifun) = czero
!    mi = abs(mval(ifun))
     mi = mval(ifun)
     do l = abs(mi), lmax1
        do ilat = 1, nlat
           do irad = llr, ulr
              orb2(irad, l, ifun) = orb2(irad, l,    ifun) &
                                  + orb1(irad, ilat, ifun) * legb1(ilat, l, mi)
           end do
        end do
     end do
!     ! Condon-Shortley factor
!     if (mval(ifun) > 0 .and. mod(mval(ifun),2) == 1) then
!        do l = abs(mi), lmax1
!           do irad = llr, ulr
!              orb2(irad, l, ifun) = - orb2(irad, l, ifun)
!           end do
!        end do
!     end if
  end do
  !$omp end parallel

end subroutine bas_ang2sph1_dyn_radfc
!######################################################################
subroutine bas_ang2sph1_fc(orb1, orb2)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad, nradfc
  use mod_sph, only : lmax1, nlat, legb1
  use mod_bas, only : mval
  use mod_ormas, only : nfcore
  use mod_const, only : czero

  implicit none
  complex(c_double_complex), intent(in) :: orb1(1:(nrad-1), 1:nlat, 1:*)
  complex(c_double_complex), intent(out) :: orb2(1:(nrad-1), 0:lmax1, 1:*)

  integer(c_long) :: ifun, irad, ilat, l, mi, llr, ulr

  !$omp parallel default(shared) private(mi, llr, ulr)
! call util_omp_disp(1, nrad - 1, llr, ulr)
  call util_omp_disp(1, nradfc, llr, ulr)

  do ifun = 1, nfcore
     orb2(llr:ulr, 0:lmax1, ifun) = czero
!    mi = abs(mval(ifun))
     mi = mval(ifun)
     do l = abs(mi), lmax1
        do ilat = 1, nlat
           do irad = llr, ulr
              orb2(irad, l, ifun) = orb2(irad, l,    ifun) &
                                  + orb1(irad, ilat, ifun) * legb1(ilat, l, mi)
           end do
        end do
     end do
!     ! Condon-Shortley factor
!     if (mval(ifun) > 0 .and. mod(mval(ifun),2) == 1) then
!        do l = abs(mi), lmax1
!           do irad = llr, ulr
!              orb2(irad, l, ifun) = - orb2(irad, l, ifun)
!           end do
!        end do
!     end if
  end do
  !$omp end parallel

end subroutine bas_ang2sph1_fc
!old!######################################################################
!oldsubroutine bas_ang2sph1_fc1(orb1, orb2)
!old
!old  use, intrinsic :: iso_c_binding
!old  use mod_rad, only : nrad, nradfc
!old  use mod_sph, only : lmax1, nlat, legb1
!old  use mod_bas, only : mval
!old  use mod_ormas, only : nfcore2, nfcore
!old  use mod_const, only : czero
!old
!old  implicit none
!old  complex(c_double_complex), intent(in) :: orb1(1:(nrad-1), 1:nlat, 1:*)
!old  complex(c_double_complex), intent(out) :: orb2(1:(nrad-1), 0:lmax1, 1:*)
!old
!old  integer(c_long) :: ifun, irad, ilat, l, mi, llr, ulr
!old
!old  !$omp parallel default(shared) private(mi, llr, ulr)
!old! call util_omp_disp(1, nrad - 1, llr, ulr)
!old  call util_omp_disp(1, nradfc, llr, ulr)
!old
!old  do ifun = nfcore2 + 1, nfcore
!old     orb2(llr:ulr, 0:lmax1, ifun) = czero
!old!    mi = abs(mval(ifun))
!old     mi = mval(ifun)
!old     do l = abs(mi), lmax1
!old        do ilat = 1, nlat
!old           do irad = llr, ulr
!old              orb2(irad, l, ifun) = orb2(irad, l,    ifun) &
!old                                  + orb1(irad, ilat, ifun) * legb1(ilat, l, mi)
!old           end do
!old        end do
!old     end do
!old!     ! Condon-Shortley factor
!old!     if (mval(ifun) > 0 .and. mod(mval(ifun),2) == 1) then
!old!        do l = abs(mi), lmax1
!old!           do irad = llr, ulr
!old!              orb2(irad, l, ifun) = - orb2(irad, l, ifun)
!old!           end do
!old!        end do
!old!     end if
!old  end do
!old  !$omp end parallel
!old
!oldend subroutine bas_ang2sph1_fc1
!old!######################################################################
!oldsubroutine bas_ang2sph1_fc2(orb1, orb2)
!old
!old  use, intrinsic :: iso_c_binding
!old  use mod_rad, only : nrad, nradfc
!old  use mod_sph, only : lmax1, nlat, legb1
!old  use mod_bas, only : mval
!old  use mod_ormas, only : nfcore2, nfcore
!old  use mod_const, only : czero
!old
!old  implicit none
!old  complex(c_double_complex), intent(in) :: orb1(1:(nrad-1), 1:nlat, 1:*)
!old  complex(c_double_complex), intent(out) :: orb2(1:(nrad-1), 0:lmax1, 1:*)
!old
!old  integer(c_long) :: ifun, irad, ilat, l, mi, llr, ulr
!old
!old  !$omp parallel default(shared) private(mi, llr, ulr)
!old! call util_omp_disp(1, nrad - 1, llr, ulr)
!old  call util_omp_disp(1, nradfc, llr, ulr)
!old
!old  do ifun = 1, nfcore2
!old     orb2(llr:ulr, 0:lmax1, ifun) = czero
!old!    mi = abs(mval(ifun))
!old     mi = mval(ifun)
!old     do l = abs(mi), lmax1
!old        do ilat = 1, nlat
!old           do irad = llr, ulr
!old              orb2(irad, l, ifun) = orb2(irad, l,    ifun) &
!old                                  + orb1(irad, ilat, ifun) * legb1(ilat, l, mi)
!old           end do
!old        end do
!old     end do
!old!     ! Condon-Shortley factor
!old!     if (mval(ifun) > 0 .and. mod(mval(ifun),2) == 1) then
!old!        do l = abs(mi), lmax1
!old!           do irad = llr, ulr
!old!              orb2(irad, l, ifun) = - orb2(irad, l, ifun)
!old!           end do
!old!        end do
!old!     end if
!old  end do
!old  !$omp end parallel
!old
!oldend subroutine bas_ang2sph1_fc2
!old!######################################################################
