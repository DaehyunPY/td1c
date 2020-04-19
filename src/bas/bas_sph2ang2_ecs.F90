!######################################################################
subroutine bas_sph2ang2_ecs(fun1, fun2)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad
  use mod_sph, only : lmax2, nlat, legf2
  use mod_bas, only : mval
  use mod_ormas, only : nfun
  use mod_const, only : czero

  implicit none
  complex(c_double_complex), intent(in) :: fun1(1:(nrad-1), 0:lmax2, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(out) :: fun2(1:(nrad-1), 1:nlat, 1:nfun, 1:nfun)

  integer(c_int) :: ifun, jfun, irad, ilat, l, mji, llr, ulr

  !$omp parallel default(shared) private(mji, llr, ulr)
  call util_omp_disp(1, nrad - 1, llr, ulr)

  do ifun = 1, nfun
     do jfun = 1, nfun !! <<<< changed
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
     end do
  end do
  !$omp end parallel

end subroutine bas_sph2ang2_ecs
!######################################################################
subroutine bas_sph2ang2_dyn_ecs(fun1, fun2)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad
  use mod_sph, only : lmax2, nlat, legf2
  use mod_bas, only : mval
  use mod_ormas, only : nfcore, nfun
  use mod_const, only : czero

  implicit none
  complex(c_double_complex), intent(in) :: fun1(1:(nrad-1), 0:lmax2, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(out) :: fun2(1:(nrad-1), 1:nlat, 1:nfun, 1:nfun)

  integer(c_int) :: ifun, jfun, irad, ilat, l, mji, llr, ulr


  integer(c_int) :: ijl, nijl, mapijl(1:3, 1:nfun*nfun*nlat)
!  real(c_double) :: xvec(1:(nrad-1), 1:2), yvec(1:(nrad-1), 1:2)
  nijl = 0
  do ifun = nfcore + 1, nfun
     do jfun = nfcore + 1, nfun !! <<<< changed
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

!DEBUG  integer(c_int) :: ipair, npair, mappair(1:2, 1:nfun*nfun)
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

end subroutine bas_sph2ang2_dyn_ecs
!######################################################################
subroutine bas_sph2ang2_fcdyn_ecs(fun1, fun2)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad, nradfc
  use mod_sph, only : lmax2, nlat, legf2
  use mod_bas, only : mval
  use mod_ormas, only : nfcore, nfun
  use mod_const, only : czero

  implicit none
  complex(c_double_complex), intent(in) :: fun1(1:(nrad-1), 0:lmax2, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(out) :: fun2(1:(nrad-1), 1:nlat, 1:nfun, 1:nfun)

  integer(c_int) :: ifun, jfun, irad, ilat, l, mji, llr, ulr

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

end subroutine bas_sph2ang2_fcdyn_ecs
!######################################################################
subroutine bas_sph2ang2_fc1dyn_ecs(fun1, fun2)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad, nradfc
  use mod_sph, only : lmax2, nlat, legf2
  use mod_bas, only : mval
  use mod_ormas, only : nfcore2, nfcore, nfun
  use mod_const, only : czero

  implicit none
  complex(c_double_complex), intent(in) :: fun1(1:(nrad-1), 0:lmax2, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(out) :: fun2(1:(nrad-1), 1:nlat, 1:nfun, 1:nfun)

  integer(c_int) :: ifun, jfun, irad, ilat, l, mji, llr, ulr

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
     end do
  end do

  do ifun = nfcore + 1, nfun
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
     end do
  end do
  !$omp end parallel

end subroutine bas_sph2ang2_fc1dyn_ecs
!######################################################################
subroutine bas_sph2ang2_fcfc_ecs(fun1, fun2)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad
  use mod_sph, only : lmax2, nlat, legf2
  use mod_bas, only : mval
  use mod_ormas, only : nfcore, nfun
  use mod_const, only : czero

  implicit none
  complex(c_double_complex), intent(in) :: fun1(1:(nrad-1), 0:lmax2, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(out) :: fun2(1:(nrad-1), 1:nlat, 1:nfun, 1:nfun)

  integer(c_int) :: ifun, jfun, irad, ilat, l, mji, llr, ulr

  !$omp parallel default(shared) private(mji, llr, ulr)
  call util_omp_disp(1, nrad - 1, llr, ulr)

  do ifun = 1, nfcore
     do jfun = 1, nfcore ! <<<<<<<< changed
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
     end do
  end do
  !$omp end parallel

end subroutine bas_sph2ang2_fcfc_ecs
!######################################################################
