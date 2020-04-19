!######################################################################
subroutine bas_sph2ang2(fun1, fun2)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad
  use mod_sph, only : lmax2, nlat, legf2
  use mod_bas, only : mval
  use mod_ormas, only : nfun
  use mod_const, only : czero

  implicit none
  complex(c_double_complex), intent(in) :: fun1(1:(nrad-1), 0:lmax2, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(out) :: fun2(1:(nrad-1), 1:nlat, 1:nfun, 1:nfun)

  integer(c_long) :: ifun, jfun, irad, ilat, l, mji, llr, ulr

  !$omp parallel default(shared) private(mji, llr, ulr)
  call util_omp_disp(1, nrad - 1, llr, ulr)

  do ifun = 1, nfun
     do jfun = ifun, nfun
        fun2(llr:ulr, 1:nlat, jfun, ifun) = czero

!       mji = abs(-mval(jfun) + mval(ifun))
        mji = -mval(jfun) + mval(ifun)
        do ilat = 1, nlat
           do l = abs(mji), lmax2
              do irad = llr, ulr
                 fun2(irad, ilat, jfun, ifun) = fun2(irad, ilat, jfun, ifun) &
                                              + fun1(irad, l,    jfun, ifun) * legf2(l, ilat, mji)
              end do
           end do
        end do
!        ! Condon-Shortley factor
!        if (-mval(jfun)+mval(ifun) > 0 .and. mod(-mval(jfun)+mval(ifun),2) == 1) then
!           do ilat = 1, nlat
!              do irad = llr, ulr
!                 fun2(irad, ilat, jfun, ifun) = - fun2(irad, ilat, jfun, ifun)
!              end do
!           end do
!        end if
     end do
  end do
  !$omp end parallel

end subroutine bas_sph2ang2
!######################################################################
subroutine bas_sph2ang2_fc(fun1, fun2)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad
  use mod_sph, only : lmax2, nlat, legf2
  use mod_bas, only : mval
  use mod_ormas, only : nfcore, nfun
  use mod_const, only : czero

  implicit none
  complex(c_double_complex), intent(in) :: fun1(1:(nrad-1), 0:lmax2, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(out) :: fun2(1:(nrad-1), 1:nlat, 1:nfun, 1:nfun)
  integer(c_long) :: ifun, irad, ilat, l, mji, llr, ulr

  !$omp parallel default(shared) private(mji, llr, ulr)
  call util_omp_disp(1, nrad - 1, llr, ulr)

  do ifun = 1, nfcore
     fun2(llr:ulr, 1:nlat, ifun, ifun) = czero

     mji = 0
     do ilat = 1, nlat
        do l = abs(mji), lmax2
           do irad = llr, ulr
              fun2(irad, ilat, ifun, ifun) = fun2(irad, ilat, ifun, ifun) &
                                           + fun1(irad, l,    ifun, ifun) * legf2(l, ilat, mji)
           end do
        end do
     end do
  end do
  !$omp end parallel

end subroutine bas_sph2ang2_fc
!######################################################################
subroutine bas_sph2ang2_fc1(fun1, fun2)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad
  use mod_sph, only : lmax2, nlat, legf2
  use mod_bas, only : mval
  use mod_ormas, only : nfcore2, nfcore, nfun
  use mod_const, only : czero

  implicit none
  complex(c_double_complex), intent(in) :: fun1(1:(nrad-1), 0:lmax2, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(out) :: fun2(1:(nrad-1), 1:nlat, 1:nfun, 1:nfun)
  integer(c_long) :: ifun, irad, ilat, l, mji, llr, ulr

  !$omp parallel default(shared) private(mji, llr, ulr)
  call util_omp_disp(1, nrad - 1, llr, ulr)

  do ifun = nfcore2 + 1, nfcore
     fun2(llr:ulr, 1:nlat, ifun, ifun) = czero

     mji = 0
     do ilat = 1, nlat
        do l = abs(mji), lmax2
           do irad = llr, ulr
              fun2(irad, ilat, ifun, ifun) = fun2(irad, ilat, ifun, ifun) &
                                           + fun1(irad, l,    ifun, ifun) * legf2(l, ilat, mji)
           end do
        end do
     end do
  end do
  !$omp end parallel

end subroutine bas_sph2ang2_fc1
!######################################################################
subroutine bas_sph2ang2_fc2(fun1, fun2)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad
  use mod_sph, only : lmax2, nlat, legf2
  use mod_bas, only : mval
  use mod_ormas, only : nfcore2, nfcore, nfun
  use mod_const, only : czero

  implicit none
  complex(c_double_complex), intent(in) :: fun1(1:(nrad-1), 0:lmax2, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(out) :: fun2(1:(nrad-1), 1:nlat, 1:nfun, 1:nfun)
  integer(c_long) :: ifun, irad, ilat, l, mji, llr, ulr

  !$omp parallel default(shared) private(mji, llr, ulr)
  call util_omp_disp(1, nrad - 1, llr, ulr)

  do ifun = 1, nfcore2
     fun2(llr:ulr, 1:nlat, ifun, ifun) = czero

     mji = 0
     do ilat = 1, nlat
        do l = abs(mji), lmax2
           do irad = llr, ulr
              fun2(irad, ilat, ifun, ifun) = fun2(irad, ilat, ifun, ifun) &
                                           + fun1(irad, l,    ifun, ifun) * legf2(l, ilat, mji)
           end do
        end do
     end do
  end do
  !$omp end parallel

end subroutine bas_sph2ang2_fc2
!######################################################################
subroutine bas_sph2ang2_fcfc(fun1, fun2)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad
  use mod_sph, only : lmax2, nlat, legf2
  use mod_bas, only : mval
  use mod_ormas, only : nfcore, nfun
  use mod_const, only : czero

  implicit none
  complex(c_double_complex), intent(in) :: fun1(1:(nrad-1), 0:lmax2, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(out) :: fun2(1:(nrad-1), 1:nlat, 1:nfun, 1:nfun)

  integer(c_long) :: ifun, jfun, irad, ilat, l, mji, llr, ulr

  !$omp parallel default(shared) private(mji, llr, ulr)
  call util_omp_disp(1, nrad - 1, llr, ulr)

  do ifun = 1, nfcore
     do jfun = ifun, nfcore
        fun2(llr:ulr, 1:nlat, jfun, ifun) = czero

!       mji = abs(-mval(jfun) + mval(ifun))
        mji = -mval(jfun) + mval(ifun)
        do ilat = 1, nlat
           do l = abs(mji), lmax2
              do irad = llr, ulr
                 fun2(irad, ilat, jfun, ifun) = fun2(irad, ilat, jfun, ifun) &
                                              + fun1(irad, l,    jfun, ifun) * legf2(l, ilat, mji)
              end do
           end do
        end do
!        ! Condon-Shortley factor
!        if (-mval(jfun)+mval(ifun) > 0 .and. mod(-mval(jfun)+mval(ifun),2) == 1) then
!           do ilat = 1, nlat
!              do irad = llr, ulr
!                 fun2(irad, ilat, jfun, ifun) = - fun2(irad, ilat, jfun, ifun)
!              end do
!           end do
!        end if
     end do
  end do
  !$omp end parallel

end subroutine bas_sph2ang2_fcfc
!######################################################################
subroutine bas_sph2ang2_fc1fc1(fun1, fun2)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad
  use mod_sph, only : lmax2, nlat, legf2
  use mod_bas, only : mval
  use mod_ormas, only : nfcore, nfcore2, nfun
  use mod_const, only : czero

  implicit none
  complex(c_double_complex), intent(in) :: fun1(1:(nrad-1), 0:lmax2, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(out) :: fun2(1:(nrad-1), 1:nlat, 1:nfun, 1:nfun)

  integer(c_long) :: ifun, jfun, irad, ilat, l, mji, llr, ulr

  !$omp parallel default(shared) private(mji, llr, ulr)
  call util_omp_disp(1, nrad - 1, llr, ulr)

  do ifun = nfcore2 + 1, nfcore
     do jfun = ifun, nfcore
        fun2(llr:ulr, 1:nlat, jfun, ifun) = czero

!       mji = abs(-mval(jfun) + mval(ifun))
        mji = -mval(jfun) + mval(ifun)
        do ilat = 1, nlat
           do l = abs(mji), lmax2
              do irad = llr, ulr
                 fun2(irad, ilat, jfun, ifun) = fun2(irad, ilat, jfun, ifun) &
                                              + fun1(irad, l,    jfun, ifun) * legf2(l, ilat, mji)
              end do
           end do
        end do
!        ! Condon-Shortley factor
!        if (-mval(jfun)+mval(ifun) > 0 .and. mod(-mval(jfun)+mval(ifun),2) == 1) then
!           do ilat = 1, nlat
!              do irad = llr, ulr
!                 fun2(irad, ilat, jfun, ifun) = - fun2(irad, ilat, jfun, ifun)
!              end do
!           end do
!        end if
     end do
  end do
  !$omp end parallel

end subroutine bas_sph2ang2_fc1fc1
!######################################################################
subroutine bas_sph2ang2_fc2fc2(fun1, fun2)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad
  use mod_sph, only : lmax2, nlat, legf2
  use mod_bas, only : mval
  use mod_ormas, only : nfcore, nfcore2, nfun
  use mod_const, only : czero

  implicit none
  complex(c_double_complex), intent(in) :: fun1(1:(nrad-1), 0:lmax2, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(out) :: fun2(1:(nrad-1), 1:nlat, 1:nfun, 1:nfun)

  integer(c_long) :: ifun, jfun, irad, ilat, l, mji, llr, ulr

  !$omp parallel default(shared) private(mji, llr, ulr)
  call util_omp_disp(1, nrad - 1, llr, ulr)

  do ifun = 1, nfcore2
     do jfun = ifun, nfcore2
        fun2(llr:ulr, 1:nlat, jfun, ifun) = czero

!       mji = abs(-mval(jfun) + mval(ifun))
        mji = -mval(jfun) + mval(ifun)
        do ilat = 1, nlat
           do l = abs(mji), lmax2
              do irad = llr, ulr
                 fun2(irad, ilat, jfun, ifun) = fun2(irad, ilat, jfun, ifun) &
                                              + fun1(irad, l,    jfun, ifun) * legf2(l, ilat, mji)
              end do
           end do
        end do
!        ! Condon-Shortley factor
!        if (-mval(jfun)+mval(ifun) > 0 .and. mod(-mval(jfun)+mval(ifun),2) == 1) then
!           do ilat = 1, nlat
!              do irad = llr, ulr
!                 fun2(irad, ilat, jfun, ifun) = - fun2(irad, ilat, jfun, ifun)
!              end do
!           end do
!        end if
     end do
  end do
  !$omp end parallel

end subroutine bas_sph2ang2_fc2fc2
!######################################################################
subroutine bas_sph2ang2_fc2fc1(fun1, fun2)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad
  use mod_sph, only : lmax2, nlat, legf2
  use mod_bas, only : mval
  use mod_ormas, only : nfcore, nfcore2, nfun
  use mod_const, only : czero

  implicit none
  complex(c_double_complex), intent(in) :: fun1(1:(nrad-1), 0:lmax2, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(out) :: fun2(1:(nrad-1), 1:nlat, 1:nfun, 1:nfun)

  integer(c_long) :: ifun, jfun, irad, ilat, l, mji, llr, ulr

  !$omp parallel default(shared) private(mji, llr, ulr)
  call util_omp_disp(1, nrad - 1, llr, ulr)

  do ifun = 1, nfcore2
     do jfun = nfcore2 + 1, nfcore
        fun2(llr:ulr, 1:nlat, jfun, ifun) = czero

!       mji = abs(-mval(jfun) + mval(ifun))
        mji = -mval(jfun) + mval(ifun)
        do ilat = 1, nlat
           do l = abs(mji), lmax2
              do irad = llr, ulr
                 fun2(irad, ilat, jfun, ifun) = fun2(irad, ilat, jfun, ifun) &
                                              + fun1(irad, l,    jfun, ifun) * legf2(l, ilat, mji)
              end do
           end do
        end do
!        ! Condon-Shortley factor
!        if (-mval(jfun)+mval(ifun) > 0 .and. mod(-mval(jfun)+mval(ifun),2) == 1) then
!           do ilat = 1, nlat
!              do irad = llr, ulr
!                 fun2(irad, ilat, jfun, ifun) = - fun2(irad, ilat, jfun, ifun)
!              end do
!           end do
!        end if
     end do
  end do
  !$omp end parallel

end subroutine bas_sph2ang2_fc2fc1
!######################################################################
subroutine bas_sph2ang2_dyn(fun1, fun2)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad
  use mod_sph, only : lmax2, nlat, legf2
  use mod_bas, only : mval
  use mod_ormas, only : nfcore, nfun
  use mod_const, only : czero

  implicit none
  complex(c_double_complex), intent(in) :: fun1(1:(nrad-1), 0:lmax2, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(out) :: fun2(1:(nrad-1), 1:nlat, 1:nfun, 1:nfun)

  integer(c_long) :: ifun, jfun, irad, ilat, l, mji, llr, ulr


  integer(c_long) :: ijl, nijl, mapijl(1:3, 1:nfun*nfun*nlat)
!  real(c_double) :: xvec(1:(nrad-1), 1:2), yvec(1:(nrad-1), 1:2)
  nijl = 0
  do ifun = nfcore + 1, nfun
     do jfun = ifun, nfun
        do ilat = 1, nlat
           nijl = nijl + 1
           mapijl(1, nijl) = ifun
           mapijl(2, nijl) = jfun
           mapijl(3, nijl) = ilat
        end do
     end do
  end do

  !$omp parallel default(shared) private(mji,ifun,jfun,ilat)
  !$omp do
  do ijl = 1, nijl
     ifun = mapijl(1, ijl)
     jfun = mapijl(2, ijl)
     ilat = mapijl(3, ijl)
     mji = -mval(jfun) + mval(ifun)
!     yvec(1:(nrad-1), 1:2) = zero
!     do l = abs(mji), lmax2
!        xvec(1:(nrad-1), 1) = dble (fun1(1:(nrad-1), l, jfun, ifun))
!        xvec(1:(nrad-1), 2) = aimag(fun1(1:(nrad-1), l, jfun, ifun))
!        call daxpy(nrad-1, legf2(l, ilat, mji), xvec(1,1), 1, yvec(1,1), 1)
!        call daxpy(nrad-1, legf2(l, ilat, mji), xvec(1,2), 1, yvec(1,2), 1)
!     end do
!     fun2(1:(nrad-1), ilat, jfun, ifun) = yvec(1:(nrad-1), 1) + yvec(1:(nrad-1), 2) * iunit
     fun2(1:(nrad-1), ilat, jfun, ifun) = czero
     do l = abs(mji), lmax2
        do irad = 1, nrad - 1
           fun2(irad, ilat, jfun, ifun) = fun2(irad, ilat, jfun, ifun) &
                                        + fun1(irad, l,    jfun, ifun) * legf2(l, ilat, mji)
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel

!DEBUG  integer(c_long) :: ipair, npair, mappair(1:2, 1:nfun*nfun)
!DEBUG  npair = 0
!DEBUG  do ifun = nfcore + 1, nfun
!DEBUG     do jfun = ifun, nfun
!DEBUG        npair = npair + 1
!DEBUG        mappair(1, npair) = ifun
!DEBUG        mappair(2, npair) = jfun
!DEBUG     end do
!DEBUG  end do
!DEBUG
!DEBUG  !$omp parallel default(shared) private(mji,ifun,jfun,ilat)
!DEBUG  !$omp do
!DEBUG  do ipair = 1, npair
!DEBUG     ifun = mappair(1, ipair)
!DEBUG     jfun = mappair(2, ipair)
!DEBUG!  do ifun = nfcore + 1, nfun
!DEBUG!     do jfun = ifun, nfun
!DEBUG        fun2(1:(nrad-1), 1:nlat, jfun, ifun) = czero
!DEBUG
!DEBUG!       mji = abs(-mval(jfun) + mval(ifun))
!DEBUG        mji = -mval(jfun) + mval(ifun)
!DEBUG        do ilat = 1, nlat
!DEBUG           do l = abs(mji), lmax2
!DEBUG              do irad = 1, nrad - 1
!DEBUG                 fun2(irad, ilat, jfun, ifun) = fun2(irad, ilat, jfun, ifun) &
!DEBUG                                              + fun1(irad, l,    jfun, ifun) * legf2(l, ilat, mji)
!DEBUG              end do
!DEBUG           end do
!DEBUG        end do
!DEBUG!        ! Condon-Shortley factor
!DEBUG!        if (-mval(jfun)+mval(ifun) > 0 .and. mod(-mval(jfun)+mval(ifun),2) == 1) then
!DEBUG!           do ilat = 1, nlat
!DEBUG!              do irad = 1, nrad - 1
!DEBUG!                 fun2(irad, ilat, jfun, ifun) = - fun2(irad, ilat, jfun, ifun)
!DEBUG!              end do
!DEBUG!           end do
!DEBUG!        end if
!DEBUG
!DEBUG!     end do
!DEBUG!  end do
!DEBUG  end do
!DEBUG  !$omp end do  
!DEBUG  !$omp end parallel

!DEBUG  !$omp parallel default(shared) private(mji, llr, ulr)
!DEBUG  call util_omp_disp(1, nrad - 1, llr, ulr)
!DEBUG
!DEBUG  do ifun = nfcore + 1, nfun
!DEBUG     do jfun = ifun, nfun
!DEBUG        fun2(llr:ulr, 1:nlat, jfun, ifun) = czero
!DEBUG
!DEBUG!       mji = abs(-mval(jfun) + mval(ifun))
!DEBUG        mji = -mval(jfun) + mval(ifun)
!DEBUG        do ilat = 1, nlat
!DEBUG           do l = abs(mji), lmax2
!DEBUG              do irad = llr, ulr
!DEBUG                 fun2(irad, ilat, jfun, ifun) = fun2(irad, ilat, jfun, ifun) &
!DEBUG                                              + fun1(irad, l,    jfun, ifun) * legf2(l, ilat, mji)
!DEBUG              end do
!DEBUG           end do
!DEBUG        end do
!DEBUG!        ! Condon-Shortley factor
!DEBUG!        if (-mval(jfun)+mval(ifun) > 0 .and. mod(-mval(jfun)+mval(ifun),2) == 1) then
!DEBUG!           do ilat = 1, nlat
!DEBUG!              do irad = llr, ulr
!DEBUG!                 fun2(irad, ilat, jfun, ifun) = - fun2(irad, ilat, jfun, ifun)
!DEBUG!              end do
!DEBUG!           end do
!DEBUG!        end if
!DEBUG     end do
!DEBUG  end do
!DEBUG  !$omp end parallel

end subroutine bas_sph2ang2_dyn
!######################################################################
subroutine bas_sph2ang2_dyn3j(fun1, fun2e, fun2o)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad
  use mod_sph, only : lmax2, nlat, legf2
  use mod_bas, only : mval
  use mod_ormas, only : nfcore, nfun
  use mod_const, only : czero

  implicit none
  complex(c_double_complex), intent(in) :: fun1(1:(nrad-1), 0:lmax2, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(out) :: fun2e(1:(nrad-1), 1:nlat, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(out) :: fun2o(1:(nrad-1), 1:nlat, 1:nfun, 1:nfun)

  integer(c_long) :: ifun, jfun, irad, ilat, l, mji, llr, ulr


  integer(c_long) :: ijl, nijl, mapijl(1:3, 1:nfun*nfun*nlat)
  nijl = 0
  do ifun = nfcore + 1, nfun
     do jfun = ifun, nfun
        do ilat = 1, nlat
           nijl = nijl + 1
           mapijl(1, nijl) = ifun
           mapijl(2, nijl) = jfun
           mapijl(3, nijl) = ilat
        end do
     end do
  end do

  !$omp parallel default(shared) private(mji,ifun,jfun,ilat)
  !$omp do
  do ijl = 1, nijl
     ifun = mapijl(1, ijl)
     jfun = mapijl(2, ijl)
     ilat = mapijl(3, ijl)
     mji = -mval(jfun) + mval(ifun)
     fun2e(1:(nrad-1), ilat, jfun, ifun) = czero
     fun2o(1:(nrad-1), ilat, jfun, ifun) = czero
     do l = abs(mji), lmax2
        if (mod(l,2) == 0) then
           do irad = 1, nrad - 1
              fun2e(irad, ilat, jfun, ifun) = fun2e(irad, ilat, jfun, ifun) &
                                            + fun1 (irad, l,    jfun, ifun) * legf2(l, ilat, mji)
           end do
        else
           do irad = 1, nrad - 1
              fun2o(irad, ilat, jfun, ifun) = fun2o(irad, ilat, jfun, ifun) &
                                            + fun1 (irad, l,    jfun, ifun) * legf2(l, ilat, mji)
           end do
        end if
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine bas_sph2ang2_dyn3j
!######################################################################
subroutine bas_sph2ang2_fcdyn(fun1, fun2)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad, nradfc
  use mod_sph, only : lmax2, nlat, legf2
  use mod_bas, only : mval
  use mod_ormas, only : nfcore, nfun
  use mod_const, only : czero

  implicit none
  complex(c_double_complex), intent(in) :: fun1(1:(nrad-1), 0:lmax2, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(out) :: fun2(1:(nrad-1), 1:nlat, 1:nfun, 1:nfun)

  integer(c_long) :: ifun, jfun, irad, ilat, l, mji, llr, ulr

  !$omp parallel default(shared) private(mji, llr, ulr)
  call util_omp_disp(1, nradfc, llr, ulr)

  do ifun = 1, nfcore
     do jfun = nfcore + 1, nfun
        fun2(llr:ulr, 1:nlat, jfun, ifun) = czero

!       mji = abs(-mval(jfun) + mval(ifun))
        mji = -mval(jfun) + mval(ifun)
        do ilat = 1, nlat
           do l = abs(mji), lmax2
              do irad = llr, ulr
                 fun2(irad, ilat, jfun, ifun) = fun2(irad, ilat, jfun, ifun) &
                                              + fun1(irad, l,    jfun, ifun) * legf2(l, ilat, mji)
              end do
           end do
        end do
!        ! Condon-Shortley factor
!        if (-mval(jfun)+mval(ifun) > 0 .and. mod(-mval(jfun)+mval(ifun),2) == 1) then
!           do ilat = 1, nlat
!              do irad = llr, ulr
!                 fun2(irad, ilat, jfun, ifun) = - fun2(irad, ilat, jfun, ifun)
!              end do
!           end do
!        end if
     end do
  end do
  !$omp end parallel

end subroutine bas_sph2ang2_fcdyn
!######################################################################
subroutine bas_sph2ang2_fc1dyn(fun1, fun2)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad, nradfc
  use mod_sph, only : lmax2, nlat, legf2
  use mod_bas, only : mval
  use mod_ormas, only : nfcore2, nfcore, nfun
  use mod_const, only : czero

  implicit none
  complex(c_double_complex), intent(in) :: fun1(1:(nrad-1), 0:lmax2, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(out) :: fun2(1:(nrad-1), 1:nlat, 1:nfun, 1:nfun)

  integer(c_long) :: ifun, jfun, irad, ilat, l, mji, llr, ulr

  !$omp parallel default(shared) private(mji, llr, ulr)
  call util_omp_disp(1, nradfc, llr, ulr)

  do ifun = nfcore2 + 1, nfcore
     do jfun = nfcore + 1, nfun
        fun2(llr:ulr, 1:nlat, jfun, ifun) = czero

!       mji = abs(-mval(jfun) + mval(ifun))
        mji = -mval(jfun) + mval(ifun)
        do ilat = 1, nlat
           do l = abs(mji), lmax2
              do irad = llr, ulr
                 fun2(irad, ilat, jfun, ifun) = fun2(irad, ilat, jfun, ifun) &
                                              + fun1(irad, l,    jfun, ifun) * legf2(l, ilat, mji)
              end do
           end do
        end do
!        ! Condon-Shortley factor
!        if (-mval(jfun)+mval(ifun) > 0 .and. mod(-mval(jfun)+mval(ifun),2) == 1) then
!           do ilat = 1, nlat
!              do irad = llr, ulr
!                 fun2(irad, ilat, jfun, ifun) = - fun2(irad, ilat, jfun, ifun)
!              end do
!           end do
!        end if
     end do
  end do
  !$omp end parallel

end subroutine bas_sph2ang2_fc1dyn
!######################################################################
subroutine bas_sph2ang2_fc1dyn3j(fun1, fun2e, fun2o)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad, nradfc
  use mod_sph, only : lmax2, nlat, legf2
  use mod_bas, only : mval
  use mod_ormas, only : nfcore2, nfcore, nfun
  use mod_const, only : czero

  implicit none
  complex(c_double_complex), intent(in) :: fun1(1:(nrad-1), 0:lmax2, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(out) :: fun2e(1:(nrad-1), 1:nlat, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(out) :: fun2o(1:(nrad-1), 1:nlat, 1:nfun, 1:nfun)

  integer(c_long) :: ifun, jfun, irad, ilat, l, mji, llr, ulr

  !$omp parallel default(shared) private(mji, llr, ulr)
  call util_omp_disp(1, nradfc, llr, ulr)

  do ifun = nfcore2 + 1, nfcore
     do jfun = nfcore + 1, nfun
        fun2e(llr:ulr, 1:nlat, jfun, ifun) = czero
        fun2o(llr:ulr, 1:nlat, jfun, ifun) = czero

!       mji = abs(-mval(jfun) + mval(ifun))
        mji = -mval(jfun) + mval(ifun)
        do ilat = 1, nlat
           do l = abs(mji), lmax2
              if (mod(l,2) == 0) then
                 do irad = llr, ulr
                    fun2e(irad, ilat, jfun, ifun) = fun2e(irad, ilat, jfun, ifun) &
                                                  + fun1 (irad, l,    jfun, ifun) * legf2(l, ilat, mji)
                 end do
              else
                 do irad = llr, ulr
                    fun2o(irad, ilat, jfun, ifun) = fun2o(irad, ilat, jfun, ifun) &
                                                  + fun1 (irad, l,    jfun, ifun) * legf2(l, ilat, mji)
                 end do
              end if
           end do
        end do
!        ! Condon-Shortley factor
!        if (-mval(jfun)+mval(ifun) > 0 .and. mod(-mval(jfun)+mval(ifun),2) == 1) then
!           do ilat = 1, nlat
!              do irad = llr, ulr
!                 fun2(irad, ilat, jfun, ifun) = - fun2(irad, ilat, jfun, ifun)
!              end do
!           end do
!        end if
     end do
  end do
  !$omp end parallel

end subroutine bas_sph2ang2_fc1dyn3j
!######################################################################
subroutine bas_sph2ang2_fc2dyn(fun1, fun2)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad, nradfc
  use mod_sph, only : lmax2, nlat, legf2
  use mod_bas, only : mval
  use mod_ormas, only : nfcore2, nfcore, nfun
  use mod_const, only : czero

  implicit none
  complex(c_double_complex), intent(in) :: fun1(1:(nrad-1), 0:lmax2, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(out) :: fun2(1:(nrad-1), 1:nlat, 1:nfun, 1:nfun)

  integer(c_long) :: ifun, jfun, irad, ilat, l, mji, llr, ulr

  !$omp parallel default(shared) private(mji, llr, ulr)
  call util_omp_disp(1, nradfc, llr, ulr)

  do ifun = 1, nfcore2
     do jfun = nfcore + 1, nfun
        fun2(llr:ulr, 1:nlat, jfun, ifun) = czero

!       mji = abs(-mval(jfun) + mval(ifun))
        mji = -mval(jfun) + mval(ifun)
        do ilat = 1, nlat
           do l = abs(mji), lmax2
              do irad = llr, ulr
                 fun2(irad, ilat, jfun, ifun) = fun2(irad, ilat, jfun, ifun) &
                                              + fun1(irad, l,    jfun, ifun) * legf2(l, ilat, mji)
              end do
           end do
        end do
!        ! Condon-Shortley factor
!        if (-mval(jfun)+mval(ifun) > 0 .and. mod(-mval(jfun)+mval(ifun),2) == 1) then
!           do ilat = 1, nlat
!              do irad = llr, ulr
!                 fun2(irad, ilat, jfun, ifun) = - fun2(irad, ilat, jfun, ifun)
!              end do
!           end do
!        end if
     end do
  end do
  !$omp end parallel

end subroutine bas_sph2ang2_fc2dyn
!######################################################################
subroutine bas_sph2ang2one(valm, fun1, fun2)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad, nradfc
  use mod_sph, only : lmax2, nlat, legf2
  use mod_const, only : czero

  implicit none
  integer(c_long), intent(in) :: valm
  complex(c_double_complex), intent(in) :: fun1(1:(nrad-1), 0:lmax2)
  complex(c_double_complex), intent(out) :: fun2(1:(nrad-1), 1:nlat)

  integer(c_long) :: irad, ilat, l, llr, ulr

  !$omp parallel default(shared) private(llr, ulr)
  call util_omp_disp(1, nradfc, llr, ulr)

  fun2(llr:ulr, 1:nlat) = czero
  do ilat = 1, nlat
     do l = abs(valm), lmax2
        do irad = llr, ulr
           fun2(irad, ilat) = fun2(irad, ilat) &
                            + fun1(irad, l) * legf2(l, ilat, valm)
        end do
     end do
  end do
!        ! Condon-Shortley factor
!        if (-mval(jfun)+mval(ifun) > 0 .and. mod(-mval(jfun)+mval(ifun),2) == 1) then
!           do ilat = 1, nlat
!              do irad = llr, ulr
!                 fun2(irad, ilat, jfun, ifun) = - fun2(irad, ilat, jfun, ifun)
!              end do
!           end do
!        end if
  !$omp end parallel

end subroutine bas_sph2ang2one
!######################################################################
