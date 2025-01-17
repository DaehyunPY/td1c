!######################################################################
subroutine bas_ang2sph2(fun1, fun2)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad
  use mod_sph, only : lmax2, nlat, legb2
  use mod_bas, only : mval
  use mod_ormas, only : nfun
  use mod_const, only : czero

  implicit none
  complex(c_double_complex), intent(in) :: fun1(1:(nrad-1), 1:nlat, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(out) :: fun2(1:(nrad-1), 0:lmax2, 1:nfun, 1:nfun)

  integer(c_int) :: ifun, jfun, irad, ilat, l, mji, llr, ulr

  !$omp parallel default(shared) private(mji, llr, ulr)
  call util_omp_disp(1, nrad - 1, llr, ulr)

  do ifun = 1, nfun
     do jfun = ifun, nfun
        fun2(llr:ulr, 0:lmax2, jfun, ifun) = czero

!       mji = abs(-mval(jfun) + mval(ifun))
        mji = -mval(jfun) + mval(ifun)
        do l = abs(mji), lmax2
           do ilat = 1, nlat
              do irad = llr, ulr
                 fun2(irad, l, jfun, ifun) = fun2(irad, l,    jfun, ifun) &
                                           + fun1(irad, ilat, jfun, ifun) * legb2(ilat, l, mji)
              end do
           end do
        end do
!        ! Condon-Shortley factor
!        if (-mval(jfun)+mval(ifun) > 0 .and. mod(-mval(jfun)+mval(ifun),2) == 1) then
!           do l = abs(mji), lmax2
!              do irad = llr, ulr
!                 fun2(irad, l, jfun, ifun) = - fun2(irad, l, jfun, ifun)
!              end do
!           end do
!        end if
     end do
  end do
  !$omp end parallel

end subroutine bas_ang2sph2
!######################################################################
subroutine bas_ang2sph2_fc(fun1, fun2)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad
  use mod_sph, only : lmax2, nlat, legb2
  use mod_bas, only : mval
  use mod_ormas, only : nfcore, nfun
  use mod_const, only : czero

  implicit none
  complex(c_double_complex), intent(in) :: fun1(1:(nrad-1), 1:nlat, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(out) :: fun2(1:(nrad-1), 0:lmax2, 1:nfun, 1:nfun)
  integer(c_int) :: ifun, irad, ilat, l, mji, llr, ulr

  !$omp parallel default(shared) private(mji, llr, ulr)
  call util_omp_disp(1, nrad - 1, llr, ulr)

  do ifun = 1, nfcore
     fun2(llr:ulr, 0:lmax2, ifun, ifun) = czero

     mji = 0
     do l = mji, lmax2
        do ilat = 1, nlat
           do irad = llr, ulr
              fun2(irad, l, ifun, ifun) = fun2(irad, l,    ifun, ifun) &
                                        + fun1(irad, ilat, ifun, ifun) * legb2(ilat, l, mji)
           end do
        end do
     end do
  end do
  !$omp end parallel

end subroutine bas_ang2sph2_fc
!######################################################################
subroutine bas_ang2sph2_fc1(fun1, fun2)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad
  use mod_sph, only : lmax2, nlat, legb2
  use mod_bas, only : mval
  use mod_ormas, only : nfcore2, nfcore, nfun
  use mod_const, only : czero

  implicit none
  complex(c_double_complex), intent(in) :: fun1(1:(nrad-1), 1:nlat, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(out) :: fun2(1:(nrad-1), 0:lmax2, 1:nfun, 1:nfun)
  integer(c_int) :: ifun, irad, ilat, l, mji, llr, ulr

  !$omp parallel default(shared) private(mji, llr, ulr)
  call util_omp_disp(1, nrad - 1, llr, ulr)

  do ifun = nfcore2 + 1, nfcore
     fun2(llr:ulr, 0:lmax2, ifun, ifun) = czero

     mji = 0
     do l = mji, lmax2
        do ilat = 1, nlat
           do irad = llr, ulr
              fun2(irad, l, ifun, ifun) = fun2(irad, l,    ifun, ifun) &
                                        + fun1(irad, ilat, ifun, ifun) * legb2(ilat, l, mji)
           end do
        end do
     end do
  end do
  !$omp end parallel

end subroutine bas_ang2sph2_fc1
!######################################################################
subroutine bas_ang2sph2_fc2(fun1, fun2)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad
  use mod_sph, only : lmax2, nlat, legb2
  use mod_bas, only : mval
  use mod_ormas, only : nfcore2, nfcore, nfun
  use mod_const, only : czero

  implicit none
  complex(c_double_complex), intent(in) :: fun1(1:(nrad-1), 1:nlat, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(out) :: fun2(1:(nrad-1), 0:lmax2, 1:nfun, 1:nfun)
  integer(c_int) :: ifun, irad, ilat, l, mji, llr, ulr

  !$omp parallel default(shared) private(mji, llr, ulr)
  call util_omp_disp(1, nrad - 1, llr, ulr)

  do ifun = 1, nfcore2
     fun2(llr:ulr, 0:lmax2, ifun, ifun) = czero

     mji = 0
     do l = mji, lmax2
        do ilat = 1, nlat
           do irad = llr, ulr
              fun2(irad, l, ifun, ifun) = fun2(irad, l,    ifun, ifun) &
                                        + fun1(irad, ilat, ifun, ifun) * legb2(ilat, l, mji)
           end do
        end do
     end do
  end do
  !$omp end parallel

end subroutine bas_ang2sph2_fc2
!######################################################################
subroutine bas_ang2sph2_fcfc(fun1, fun2)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad
  use mod_sph, only : lmax2, nlat, legb2
  use mod_bas, only : mval
  use mod_ormas, only : nfcore, nfun
  use mod_const, only : czero

  implicit none
  complex(c_double_complex), intent(in) :: fun1(1:(nrad-1), 1:nlat, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(out) :: fun2(1:(nrad-1), 0:lmax2, 1:nfun, 1:nfun)
  integer(c_int) :: ifun, jfun, irad, ilat, l, mji, llr, ulr

  !$omp parallel default(shared) private(mji, llr, ulr)
  call util_omp_disp(1, nrad - 1, llr, ulr)

  do ifun = 1, nfcore
     do jfun = ifun, nfcore
        fun2(llr:ulr, 0:lmax2, jfun, ifun) = czero

!       mji = abs(-mval(jfun) + mval(ifun))
        mji = -mval(jfun) + mval(ifun)
        do l = abs(mji), lmax2
           do ilat = 1, nlat
              do irad = llr, ulr
                 fun2(irad, l, jfun, ifun) = fun2(irad, l,    jfun, ifun) &
                                           + fun1(irad, ilat, jfun, ifun) * legb2(ilat, l, mji)
              end do
           end do
        end do
!        ! Condon-Shortley factor
!        if (-mval(jfun)+mval(ifun) > 0 .and. mod(-mval(jfun)+mval(ifun),2) == 1) then
!           do l = abs(mji), lmax2
!              do irad = llr, ulr
!                 fun2(irad, l, jfun, ifun) = - fun2(irad, l, jfun, ifun)
!              end do
!           end do
!        end if
     end do
  end do
  !$omp end parallel

end subroutine bas_ang2sph2_fcfc
!######################################################################
subroutine bas_ang2sph2_fc1fc1(fun1, fun2)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad
  use mod_sph, only : lmax2, nlat, legb2
  use mod_bas, only : mval
  use mod_ormas, only : nfcore, nfcore2, nfun
  use mod_const, only : czero

  implicit none
  complex(c_double_complex), intent(in) :: fun1(1:(nrad-1), 1:nlat, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(out) :: fun2(1:(nrad-1), 0:lmax2, 1:nfun, 1:nfun)
  integer(c_int) :: ifun, jfun, irad, ilat, l, mji, llr, ulr

  !$omp parallel default(shared) private(mji, llr, ulr)
  call util_omp_disp(1, nrad - 1, llr, ulr)

  do ifun = nfcore2 + 1, nfcore
     do jfun = ifun, nfcore
        fun2(llr:ulr, 0:lmax2, jfun, ifun) = czero

!       mji = abs(-mval(jfun) + mval(ifun))
        mji = -mval(jfun) + mval(ifun)
        do l = abs(mji), lmax2
           do ilat = 1, nlat
              do irad = llr, ulr
                 fun2(irad, l, jfun, ifun) = fun2(irad, l,    jfun, ifun) &
                                           + fun1(irad, ilat, jfun, ifun) * legb2(ilat, l, mji)
              end do
           end do
        end do
!        ! Condon-Shortley factor
!        if (-mval(jfun)+mval(ifun) > 0 .and. mod(-mval(jfun)+mval(ifun),2) == 1) then
!           do l = abs(mji), lmax2
!              do irad = llr, ulr
!                 fun2(irad, l, jfun, ifun) = - fun2(irad, l, jfun, ifun)
!              end do
!           end do
!        end if
     end do
  end do
  !$omp end parallel

end subroutine bas_ang2sph2_fc1fc1
!######################################################################
subroutine bas_ang2sph2_fc2fc2(fun1, fun2)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad
  use mod_sph, only : lmax2, nlat, legb2
  use mod_bas, only : mval
  use mod_ormas, only : nfcore, nfcore2, nfun
  use mod_const, only : czero

  implicit none
  complex(c_double_complex), intent(in) :: fun1(1:(nrad-1), 1:nlat, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(out) :: fun2(1:(nrad-1), 0:lmax2, 1:nfun, 1:nfun)
  integer(c_int) :: ifun, jfun, irad, ilat, l, mji, llr, ulr

  !$omp parallel default(shared) private(mji, llr, ulr)
  call util_omp_disp(1, nrad - 1, llr, ulr)

  do ifun = 1, nfcore2
     do jfun = ifun, nfcore2
        fun2(llr:ulr, 0:lmax2, jfun, ifun) = czero

!       mji = abs(-mval(jfun) + mval(ifun))
        mji = -mval(jfun) + mval(ifun)
        do l = abs(mji), lmax2
           do ilat = 1, nlat
              do irad = llr, ulr
                 fun2(irad, l, jfun, ifun) = fun2(irad, l,    jfun, ifun) &
                                           + fun1(irad, ilat, jfun, ifun) * legb2(ilat, l, mji)
              end do
           end do
        end do
!        ! Condon-Shortley factor
!        if (-mval(jfun)+mval(ifun) > 0 .and. mod(-mval(jfun)+mval(ifun),2) == 1) then
!           do l = abs(mji), lmax2
!              do irad = llr, ulr
!                 fun2(irad, l, jfun, ifun) = - fun2(irad, l, jfun, ifun)
!              end do
!           end do
!        end if
     end do
  end do
  !$omp end parallel

end subroutine bas_ang2sph2_fc2fc2
!######################################################################
subroutine bas_ang2sph2_fc2fc1(fun1, fun2)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad
  use mod_sph, only : lmax2, nlat, legb2
  use mod_bas, only : mval
  use mod_ormas, only : nfcore, nfcore2, nfun
  use mod_const, only : czero

  implicit none
  complex(c_double_complex), intent(in) :: fun1(1:(nrad-1), 1:nlat, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(out) :: fun2(1:(nrad-1), 0:lmax2, 1:nfun, 1:nfun)
  integer(c_int) :: ifun, jfun, irad, ilat, l, mji, llr, ulr

  !$omp parallel default(shared) private(mji, llr, ulr)
  call util_omp_disp(1, nrad - 1, llr, ulr)

  do ifun = 1, nfcore2
     do jfun = nfcore2 + 1, nfcore
        fun2(llr:ulr, 0:lmax2, jfun, ifun) = czero

!       mji = abs(-mval(jfun) + mval(ifun))
        mji = -mval(jfun) + mval(ifun)
        do l = abs(mji), lmax2
           do ilat = 1, nlat
              do irad = llr, ulr
                 fun2(irad, l, jfun, ifun) = fun2(irad, l,    jfun, ifun) &
                                           + fun1(irad, ilat, jfun, ifun) * legb2(ilat, l, mji)
              end do
           end do
        end do
!        ! Condon-Shortley factor
!        if (-mval(jfun)+mval(ifun) > 0 .and. mod(-mval(jfun)+mval(ifun),2) == 1) then
!           do l = abs(mji), lmax2
!              do irad = llr, ulr
!                 fun2(irad, l, jfun, ifun) = - fun2(irad, l, jfun, ifun)
!              end do
!           end do
!        end if
     end do
  end do
  !$omp end parallel

end subroutine bas_ang2sph2_fc2fc1
!######################################################################
subroutine bas_ang2sph2_dyn(fun1, fun2)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad
  use mod_sph, only : lmax2, nlat, legb2
  use mod_bas, only : mval
  use mod_ormas, only : nfcore, nfun
  use mod_const, only : zero, czero, iunit

  implicit none
  complex(c_double_complex), intent(in) :: fun1(1:(nrad-1), 1:nlat, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(out) :: fun2(1:(nrad-1), 0:lmax2, 1:nfun, 1:nfun)
  integer(c_int) :: ifun, jfun, irad, ilat, l, mji, llr, ulr
  integer(c_int) :: ijl, nijl, mapijl(1:3, 1:nfun*nfun*(lmax2+1))
!  real(c_double) :: xvec(1:(nrad-1), 1:2, 0:lmax2), yvec(1:(nrad-1), 1:2, 0:lmax2)

  nijl = 0
  do ifun = nfcore + 1, nfun
     do jfun = ifun, nfun
        mji = -mval(jfun) + mval(ifun)
        do l = abs(mji), lmax2
           nijl = nijl + 1
           mapijl(1, nijl) = ifun
           mapijl(2, nijl) = jfun
           mapijl(3, nijl) = l
        end do
     end do
  end do

  !$omp parallel default(shared) private(mji,ifun,jfun,l)
  !$omp do
  do ijl = 1, nijl
     ifun = mapijl(1, ijl)
     jfun = mapijl(2, ijl)
     l    = mapijl(3, ijl)
     mji = -mval(jfun) + mval(ifun)
!dame     yvec(1:(nrad-1), 1:2, l) = zero
!dame     do ilat = 1, nlat
!dame        xvec(1:(nrad-1), 1, l) = dble (fun1(1:(nrad-1), ilat, jfun, ifun))
!dame        xvec(1:(nrad-1), 2, l) = aimag(fun1(1:(nrad-1), ilat, jfun, ifun))
!dame        call daxpy(nrad-1, legb2(ilat, l, mji), xvec(1,1,l), 1, yvec(1,1,l), 1)
!dame        call daxpy(nrad-1, legb2(ilat, l, mji), xvec(1,2,l), 1, yvec(1,2,l), 1)
!dame     end do
!dame     fun2(1:(nrad-1), l, jfun, ifun) = yvec(1:(nrad-1), 1, l) + yvec(1:(nrad-1), 2, l) * iunit

!dame     do irad = 1, nrad - 1
!dame        yvec(irad, 1, l) = zero
!dame        yvec(irad, 2, l) = zero
!dame     end do
!dame     do ilat = 1, nlat
!dame        do irad = 1, nrad - 1
!dame           yvec(irad, 1, l) = yvec(irad, 1, l) + dble (fun1(irad, ilat, jfun, ifun)) * legb2(ilat, l, mji)
!dame           yvec(irad, 2, l) = yvec(irad, 2, l) + aimag(fun1(irad, ilat, jfun, ifun)) * legb2(ilat, l, mji)
!dame        end do
!dame     end do
!dame     do irad = 1, nrad - 1
!dame        fun2(irad, l, jfun, ifun) = yvec(irad, 1, l) + yvec(irad, 2, l) * iunit
!dame     end do

!dame     fun2(1:(nrad-1), l, jfun, ifun) = czero
!dame     do ilat = 1, nlat
!dame        do irad = 1, nrad - 1
!dame           fun2(irad, l, jfun, ifun) = fun2(irad, l,    jfun, ifun) &
!dame                                     +  dble(fun1(irad, ilat, jfun, ifun)) * legb2(ilat, l, mji) &
!dame                                     + aimag(fun1(irad, ilat, jfun, ifun)) * legb2(ilat, l, mji) * iunit
!dame        end do
!dame     end do

     fun2(1:(nrad-1), l, jfun, ifun) = czero
     do ilat = 1, nlat
        do irad = 1, nrad - 1
           fun2(irad, l, jfun, ifun) = fun2(irad, l,    jfun, ifun) &
                                     + fun1(irad, ilat, jfun, ifun) * legb2(ilat, l, mji)
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel

!DEBUG2  integer(c_int) :: ipair, npair, mappair(1:2, 1:nfun*nfun)
!DEBUG2  npair = 0
!DEBUG2  do ifun = nfcore + 1, nfun
!DEBUG2     do jfun = ifun, nfun
!DEBUG2        npair = npair + 1
!DEBUG2        mappair(1, npair) = ifun
!DEBUG2        mappair(2, npair) = jfun
!DEBUG2     end do
!DEBUG2  end do
!DEBUG2
!DEBUG2  !$omp parallel default(shared) private(mji,ifun,jfun,l)
!DEBUG2  !$omp do
!DEBUG2  do ipair = 1, npair
!DEBUG2     ifun = mappair(1, ipair)
!DEBUG2     jfun = mappair(2, ipair)
!DEBUG2!  do ifun = nfcore + 1, nfun
!DEBUG2!     do jfun = ifun, nfun
!DEBUG2        fun2(1:(nrad-1), 0:lmax2, jfun, ifun) = czero
!DEBUG2
!DEBUG2!       mji = abs(-mval(jfun) + mval(ifun))
!DEBUG2        mji = -mval(jfun) + mval(ifun)
!DEBUG2        do l = abs(mji), lmax2
!DEBUG2           do ilat = 1, nlat
!DEBUG2              do irad = 1, nrad - 1
!DEBUG2                 fun2(irad, l, jfun, ifun) = fun2(irad, l,    jfun, ifun) &
!DEBUG2                                           + fun1(irad, ilat, jfun, ifun) * legb2(ilat, l, mji)
!DEBUG2              end do
!DEBUG2           end do
!DEBUG2        end do
!DEBUG2!        ! Condon-Shortley factor
!DEBUG2!        if (-mval(jfun)+mval(ifun) > 0 .and. mod(-mval(jfun)+mval(ifun),2) == 1) then
!DEBUG2!           do l = abs(mji), lmax2
!DEBUG2!              do irad = 1, nrad - 1
!DEBUG2!                 fun2(irad, l, jfun, ifun) = - fun2(irad, l, jfun, ifun)
!DEBUG2!              end do
!DEBUG2!           end do
!DEBUG2!        end if
!DEBUG2
!DEBUG2!     end do
!DEBUG2!  end do
!DEBUG2  end do
!DEBUG2  !$omp end do
!DEBUG2  !$omp end parallel

!DEBUG  !$omp parallel default(shared) private(mji, llr, ulr)
!DEBUG  call util_omp_disp(1, nrad - 1, llr, ulr)
!DEBUG
!DEBUG  do ifun = nfcore + 1, nfun
!DEBUG     do jfun = ifun, nfun
!DEBUG        fun2(llr:ulr, 0:lmax2, jfun, ifun) = czero
!DEBUG
!DEBUG!       mji = abs(-mval(jfun) + mval(ifun))
!DEBUG        mji = -mval(jfun) + mval(ifun)
!DEBUG        do l = abs(mji), lmax2
!DEBUG           do ilat = 1, nlat
!DEBUG              do irad = llr, ulr
!DEBUG                 fun2(irad, l, jfun, ifun) = fun2(irad, l,    jfun, ifun) &
!DEBUG                                           + fun1(irad, ilat, jfun, ifun) * legb2(ilat, l, mji)
!DEBUG              end do
!DEBUG           end do
!DEBUG        end do
!DEBUG!        ! Condon-Shortley factor
!DEBUG!        if (-mval(jfun)+mval(ifun) > 0 .and. mod(-mval(jfun)+mval(ifun),2) == 1) then
!DEBUG!           do l = abs(mji), lmax2
!DEBUG!              do irad = llr, ulr
!DEBUG!                 fun2(irad, l, jfun, ifun) = - fun2(irad, l, jfun, ifun)
!DEBUG!              end do
!DEBUG!           end do
!DEBUG!        end if
!DEBUG     end do
!DEBUG  end do
!DEBUG  !$omp end parallel

end subroutine bas_ang2sph2_dyn
!######################################################################
subroutine bas_ang2sph2_dyn3j(fun1e, fun1o, fun2)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad
  use mod_sph, only : lmax2, nlat, legb2
  use mod_bas, only : mval
  use mod_ormas, only : nfcore, nfun
  use mod_const, only : zero, czero, iunit

  implicit none
  complex(c_double_complex), intent(in) :: fun1e(1:(nrad-1), 1:nlat, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(in) :: fun1o(1:(nrad-1), 1:nlat, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(out) :: fun2(1:(nrad-1), 0:lmax2, 1:nfun, 1:nfun)
  integer(c_int) :: ifun, jfun, irad, ilat, l, mji, llr, ulr
  integer(c_int) :: ijl, nijl, mapijl(1:3, 1:nfun*nfun*(lmax2+1))
!  real(c_double) :: xvec(1:(nrad-1), 1:2, 0:lmax2), yvec(1:(nrad-1), 1:2, 0:lmax2)

  nijl = 0
  do ifun = nfcore + 1, nfun
     do jfun = ifun, nfun
        mji = -mval(jfun) + mval(ifun)
        do l = abs(mji), lmax2
           nijl = nijl + 1
           mapijl(1, nijl) = ifun
           mapijl(2, nijl) = jfun
           mapijl(3, nijl) = l
        end do
     end do
  end do

  !$omp parallel default(shared) private(mji,ifun,jfun,l)
  !$omp do
  do ijl = 1, nijl
     ifun = mapijl(1, ijl)
     jfun = mapijl(2, ijl)
     l    = mapijl(3, ijl)
     mji = -mval(jfun) + mval(ifun)

     fun2(1:(nrad-1), l, jfun, ifun) = czero
     if (mod(l,2) == 0) then
        do ilat = 1, nlat
           do irad = 1, nrad - 1
              fun2(irad, l, jfun, ifun) = fun2 (irad, l,    jfun, ifun) &
                                        + fun1e(irad, ilat, jfun, ifun) * legb2(ilat, l, mji)
           end do
        end do
     else
        do ilat = 1, nlat
           do irad = 1, nrad - 1
              fun2(irad, l, jfun, ifun) = fun2 (irad, l,    jfun, ifun) &
                                        + fun1o(irad, ilat, jfun, ifun) * legb2(ilat, l, mji)
           end do
        end do
     end if
  end do
  !$omp end do
  !$omp end parallel

end subroutine bas_ang2sph2_dyn3j
!######################################################################
subroutine bas_ang2sph2_fcdyn(fun1, fun2)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad, nradfc
  use mod_sph, only : lmax2, nlat, legb2
  use mod_bas, only : mval
  use mod_ormas, only : nfcore, nfun
  use mod_const, only : czero

  implicit none
  complex(c_double_complex), intent(in) :: fun1(1:(nrad-1), 1:nlat, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(out) :: fun2(1:(nrad-1), 0:lmax2, 1:nfun, 1:nfun)

  integer(c_int) :: ifun, jfun, irad, ilat, l, mji, llr, ulr

  !$omp parallel default(shared) private(mji, llr, ulr)
  call util_omp_disp(1, nradfc, llr, ulr)

  do ifun = 1, nfcore
     do jfun = nfcore + 1, nfun
        fun2(llr:ulr, 0:lmax2, jfun, ifun) = czero

!       mji = abs(-mval(jfun) + mval(ifun))
        mji = -mval(jfun) + mval(ifun)
        do l = abs(mji), lmax2
           do ilat = 1, nlat
              do irad = llr, ulr
                 fun2(irad, l, jfun, ifun) = fun2(irad, l,    jfun, ifun) &
                                           + fun1(irad, ilat, jfun, ifun) * legb2(ilat, l, mji)
              end do
           end do
        end do
!        ! Condon-Shortley factor
!        if (-mval(jfun)+mval(ifun) > 0 .and. mod(-mval(jfun)+mval(ifun),2) == 1) then
!           do l = abs(mji), lmax2
!              do irad = llr, ulr
!                 fun2(irad, l, jfun, ifun) = - fun2(irad, l, jfun, ifun)
!              end do
!           end do
!        end if
     end do
  end do
  !$omp end parallel

end subroutine bas_ang2sph2_fcdyn
!######################################################################
subroutine bas_ang2sph2_fc1dyn(fun1, fun2)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad, nradfc
  use mod_sph, only : lmax2, nlat, legb2
  use mod_bas, only : mval
  use mod_ormas, only : nfcore2, nfcore, nfun
  use mod_const, only : czero

  implicit none
  complex(c_double_complex), intent(in) :: fun1(1:(nrad-1), 1:nlat, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(out) :: fun2(1:(nrad-1), 0:lmax2, 1:nfun, 1:nfun)

  integer(c_int) :: ifun, jfun, irad, ilat, l, mji, llr, ulr

  !$omp parallel default(shared) private(mji, llr, ulr)
  call util_omp_disp(1, nradfc, llr, ulr)

  do ifun = nfcore2 + 1, nfcore
     do jfun = nfcore + 1, nfun
        fun2(llr:ulr, 0:lmax2, jfun, ifun) = czero

!       mji = abs(-mval(jfun) + mval(ifun))
        mji = -mval(jfun) + mval(ifun)
        do l = abs(mji), lmax2
           do ilat = 1, nlat
              do irad = llr, ulr
                 fun2(irad, l, jfun, ifun) = fun2(irad, l,    jfun, ifun) &
                                           + fun1(irad, ilat, jfun, ifun) * legb2(ilat, l, mji)
              end do
           end do
        end do
!        ! Condon-Shortley factor
!        if (-mval(jfun)+mval(ifun) > 0 .and. mod(-mval(jfun)+mval(ifun),2) == 1) then
!           do l = abs(mji), lmax2
!              do irad = llr, ulr
!                 fun2(irad, l, jfun, ifun) = - fun2(irad, l, jfun, ifun)
!              end do
!           end do
!        end if
     end do
  end do
  !$omp end parallel

end subroutine bas_ang2sph2_fc1dyn
!######################################################################
subroutine bas_ang2sph2_fc1dyn3j(fun1e, fun1o, fun2)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad, nradfc
  use mod_sph, only : lmax2, nlat, legb2
  use mod_bas, only : mval
  use mod_ormas, only : nfcore2, nfcore, nfun
  use mod_const, only : czero

  implicit none
  complex(c_double_complex), intent(in) :: fun1e(1:(nrad-1), 1:nlat, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(in) :: fun1o(1:(nrad-1), 1:nlat, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(out) :: fun2(1:(nrad-1), 0:lmax2, 1:nfun, 1:nfun)

  integer(c_int) :: ifun, jfun, irad, ilat, l, mji, llr, ulr

  !$omp parallel default(shared) private(mji, llr, ulr)
  call util_omp_disp(1, nradfc, llr, ulr)
  do ifun = nfcore2 + 1, nfcore
     do jfun = nfcore + 1, nfun
        fun2(llr:ulr, 0:lmax2, jfun, ifun) = czero
        mji = -mval(jfun) + mval(ifun)
        do l = abs(mji), lmax2
           if (mod(l,2) == 0) then
              do ilat = 1, nlat
                 do irad = llr, ulr
                    fun2(irad, l, jfun, ifun) = fun2 (irad, l,    jfun, ifun) &
                                              + fun1e(irad, ilat, jfun, ifun) * legb2(ilat, l, mji)
                 end do
              end do
           else
              do ilat = 1, nlat
                 do irad = llr, ulr
                    fun2(irad, l, jfun, ifun) = fun2 (irad, l,    jfun, ifun) &
                                              + fun1o(irad, ilat, jfun, ifun) * legb2(ilat, l, mji)
                 end do
              end do
           end if
        end do
     end do
  end do
  !$omp end parallel

end subroutine bas_ang2sph2_fc1dyn3j
!######################################################################
subroutine bas_ang2sph2_fc2dyn(fun1, fun2)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad, nradfc
  use mod_sph, only : lmax2, nlat, legb2
  use mod_bas, only : mval
  use mod_ormas, only : nfcore2, nfcore, nfun
  use mod_const, only : czero

  implicit none
  complex(c_double_complex), intent(in) :: fun1(1:(nrad-1), 1:nlat, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(out) :: fun2(1:(nrad-1), 0:lmax2, 1:nfun, 1:nfun)

  integer(c_int) :: ifun, jfun, irad, ilat, l, mji, llr, ulr

  !$omp parallel default(shared) private(mji, llr, ulr)
  call util_omp_disp(1, nradfc, llr, ulr)

  do ifun = 1, nfcore2
     do jfun = nfcore + 1, nfun
        fun2(llr:ulr, 0:lmax2, jfun, ifun) = czero

!       mji = abs(-mval(jfun) + mval(ifun))
        mji = -mval(jfun) + mval(ifun)
        do l = abs(mji), lmax2
           do ilat = 1, nlat
              do irad = llr, ulr
                 fun2(irad, l, jfun, ifun) = fun2(irad, l,    jfun, ifun) &
                                           + fun1(irad, ilat, jfun, ifun) * legb2(ilat, l, mji)
              end do
           end do
        end do
!        ! Condon-Shortley factor
!        if (-mval(jfun)+mval(ifun) > 0 .and. mod(-mval(jfun)+mval(ifun),2) == 1) then
!           do l = abs(mji), lmax2
!              do irad = llr, ulr
!                 fun2(irad, l, jfun, ifun) = - fun2(irad, l, jfun, ifun)
!              end do
!           end do
!        end if
     end do
  end do
  !$omp end parallel

end subroutine bas_ang2sph2_fc2dyn
!######################################################################
subroutine bas_ang2sph2one(valm, fun1, fun2)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad
  use mod_const, only : czero
  use mod_sph, only : lmax2, nlat, legb2

  implicit none
  integer(c_int), intent(in) :: valm
  complex(c_double_complex), intent(in) :: fun1(1:(nrad-1), 1:nlat)
  complex(c_double_complex), intent(out) :: fun2(1:(nrad-1), 0:lmax2)

  integer(c_int) :: irad, ilat, l, llr, ulr

  !$omp parallel default(shared) private(llr, ulr)
  call util_omp_disp(1, nrad - 1, llr, ulr)

  fun2(llr:ulr, 0:lmax2) = czero
  do l = abs(valm), lmax2
     do ilat = 1, nlat
        do irad = llr, ulr
           fun2(irad, l) = fun2(irad, l) &
                         + fun1(irad, ilat) * legb2(ilat, l, valm)
        end do
     end do
  end do
!        ! Condon-Shortley factor
!        if (-mval(jfun)+mval(ifun) > 0 .and. mod(-mval(jfun)+mval(ifun),2) == 1) then
!           do l = abs(valm), lmax2
!              do irad = llr, ulr
!                 fun2(irad, l, jfun, ifun) = - fun2(irad, l, jfun, ifun)
!              end do
!           end do
!        end if
  !$omp end parallel

end subroutine bas_ang2sph2one
!######################################################################
