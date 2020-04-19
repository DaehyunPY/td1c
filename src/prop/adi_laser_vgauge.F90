!######################################################################
subroutine adi_laser_vgauge(icomp, istag, dtime, lfield, wfn)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad
  use mod_sph, only : lmax1
  use mod_ormas, only : nfun

  implicit none
  integer(c_int), intent(in) :: icomp, istag
  real(c_double), intent(in) :: dtime, lfield(1:3)
  complex(c_double_complex), intent(inout) :: wfn(1:(nrad-1), 0:lmax1, 1:nfun)

  integer(c_int) :: eve_odd
  integer(c_int) :: lmaxh, lllo2, ullo2, lo2, lval

  lmaxh = lmax1 / 2
  if (istag == 0) then
     do eve_odd = 0, 1, +1
        !$omp parallel default(shared) private(lllo2, ullo2, lval)
        !###########################
        call util_omp_disp(0, lmaxh, lllo2, ullo2)
        do lo2 = lllo2, ullo2
           lval = eve_odd + 2 * lo2
           if (lval >= lmax1) cycle
           call adi_laser_vgauge_angp(icomp, dtime, lfield, wfn, lval)
           call adi_laser_vgauge_mixp(icomp, dtime, lfield, wfn, lval)
        end do
        !###########################
        !$omp end parallel
     end do
  else
     do eve_odd = 1, 0, -1
        !$omp parallel default(shared) private(lllo2, ullo2, lval)
        !###########################
        call util_omp_disp(0, lmaxh, lllo2, ullo2)
        do lo2 = lllo2, ullo2
           lval = eve_odd + 2 * lo2
           if (lval >= lmax1) cycle
           call adi_laser_vgauge_mixp(icomp, dtime, lfield, wfn, lval)
           call adi_laser_vgauge_angp(icomp, dtime, lfield, wfn, lval)
        end do
        !###########################
        !$omp end parallel
     end do
  end if

end subroutine adi_laser_vgauge
!######################################################################
subroutine adi_laser_vgauge_angp(icomp, dtime, lfield, wfn, lval)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad, xrad
  use mod_sph, only : lmax1
  use mod_bas, only : mval
  use mod_ormas, only : nfcore, nfun
  use mod_const, only : one, two, three, half, runit, iunit

  implicit none
  integer(c_int), intent(in) :: icomp, lval
  real(c_double), intent(in) :: dtime, lfield(1:3)
  complex(c_double_complex), intent(inout) :: wfn(1:(nrad-1), 0:lmax1, 1:nfun)

  integer(c_int) :: ifun, irad
  complex(c_double_complex) :: facd, facl, beta, ombd, opbd, gammp, gammm, tmp
  complex(c_double_complex) :: wfn0, wfn1, yp, ym

  facd = (lval + one) / sqrt((two * lval + one) * (two * lval + three)) 
  if (icomp == 1) then
     facl = iunit * lfield(3) * dtime * half * facd
  else
     facl = runit * lfield(3) * dtime * half * facd
  end if

  do ifun = nfcore + 1, nfun
     if (lval < abs(mval(ifun))) cycle
     beta = facl * sqrt((lval + one) ** two - mval(ifun) ** two)
     do irad = 1, nrad - 1
        tmp = beta / xrad(irad)
        ombd = one - tmp
        opbd = one + tmp
        gammp = ombd / opbd
        gammm = opbd / ombd

        wfn0 = wfn(irad, lval,     ifun)
        wfn1 = wfn(irad, lval + 1, ifun)
        yp = (wfn0 - iunit * wfn1) * gammp
        ym = (wfn0 + iunit * wfn1) * gammm
        wfn(irad, lval,     ifun) = half * (yp + ym)
        wfn(irad, lval + 1, ifun) = half * (yp - ym) * iunit
     end do
  end do

end subroutine adi_laser_vgauge_angp
!######################################################################
subroutine adi_laser_vgauge_mixp(icomp, dtime, lfield, wfn, lval)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad, ndvr, radp
  use mod_sph, only : lmax1, mmax1
  use mod_bas, only : mval
  use mod_ormas, only : nfcore, nfun
  use mod_const, only : one, two, three, half, czero, runit, iunit

  implicit none
  integer(c_int), intent(in) :: icomp, lval
  real(c_double), intent(in) :: dtime, lfield(1:3)
  complex(c_double_complex), intent(inout) :: wfn(1:(nrad-1), 0:lmax1, 1:nfun)

  integer(c_int) :: info1, info2
  integer(c_int) :: dim, ld1, ld2, ib1, ib2, jll, jul
  integer(c_int) :: ifun, mlim, m, irad, jrad, ij
  complex(c_double_complex) :: facd, facl, beta, tmp

  complex(c_double_complex), allocatable :: pwfn(:), hpwfn(:)
  complex(c_double_complex), allocatable :: mwfn(:), hmwfn(:)
  integer(c_int), allocatable :: tpiv1(:,:), tpiv2(:,:)
  complex(c_double_complex), allocatable :: tnum1(:,:,:), tnum2(:,:,:)
  complex(c_double_complex), allocatable :: tden1(:,:,:), tden2(:,:,:)

  facd = runit / sqrt((two * lval + one) * (two * lval + three))
  if (icomp == 1) then
     facl = lfield(3) * dtime * half * facd
  else
     facl = lfield(3) * dtime * half * facd
  end if

  dim = nrad - 1
  ld1 = 2 * ndvr + 1
  ld2 = 3 * ndvr + 1
  mlim = min(lval, mmax1)

  allocate(pwfn(1:dim))
  allocate(mwfn(1:dim))
  allocate(hpwfn(1:dim))
  allocate(hmwfn(1:dim))
  allocate(tpiv1(1:dim, -mlim:mlim))
  allocate(tnum1(1:ld1, 1:dim, -mlim:mlim))
  allocate(tden1(1:ld2, 1:dim, -mlim:mlim))
  allocate(tpiv2(1:dim, -mlim:mlim))
  allocate(tnum2(1:ld1, 1:dim, -mlim:mlim))
  allocate(tden2(1:ld2, 1:dim, -mlim:mlim))

  tpiv1(1:dim, -mlim:mlim) = 0
  tnum1(1:ld1, 1:dim, -mlim:mlim) = czero
  tden1(1:ld2, 1:dim, -mlim:mlim) = czero
  tpiv2(1:dim, -mlim:mlim) = 0
  tnum2(1:ld1, 1:dim, -mlim:mlim) = czero
  tden2(1:ld2, 1:dim, -mlim:mlim) = czero
  do m = -mlim, mlim
     beta = facl * sqrt((lval + one) ** two - m ** two)
     do irad = 1, dim
        ib1 = 1 + ndvr
        ib2 = 1 + ndvr * 2
        tnum1(ib1, irad, m) = runit
        tden1(ib2, irad, m) = runit
        tnum2(ib1, irad, m) = runit
        tden2(ib2, irad, m) = runit
        jll = max(1, irad - ndvr)
        jul = min(dim, irad + ndvr)
        do jrad = jll, jul
           ib1 =     ndvr + 1 + irad - jrad
           ib2 = 2 * ndvr + 1 + irad - jrad
           ij = (2 * ndvr + 1) * irad + (jrad - (irad - ndvr)) + 1 - ndvr
!
!old
!          tmp = + nabla(ij) * beta
!new
           tmp = - radp(ib1, jrad) * beta
!
           tnum1(ib1, jrad, m) = tnum1(ib1, jrad, m) - tmp
           tden1(ib2, jrad, m) = tden1(ib2, jrad, m) + tmp
           tnum2(ib1, jrad, m) = tnum2(ib1, jrad, m) + tmp
           tden2(ib2, jrad, m) = tden2(ib2, jrad, m) - tmp
        end do
     end do

     ! LU factorization
     call zgbtrf(dim, dim, ndvr, ndvr, tden1(1,1,m), ld2, tpiv1(1,m), info1)
     call zgbtrf(dim, dim, ndvr, ndvr, tden2(1,1,m), ld2, tpiv2(1,m), info2)
  end do

  do ifun = nfcore + 1, nfun
     m = mval(ifun)
     if (lval < abs(m)) cycle

     pwfn(1:dim) = czero
     mwfn(1:dim) = czero
     hpwfn(1:dim) = czero
     hmwfn(1:dim) = czero

     pwfn(1:dim) = wfn(1:dim, lval, ifun) + wfn(1:dim, lval + 1, ifun)
     mwfn(1:dim) = wfn(1:dim, lval, ifun) - wfn(1:dim, lval + 1, ifun)
!
! ZGBMM (instead of ZGBMV) and
! ZGBTRS with multiple vectors
! may be used here
!
     call ZGBMV('N', dim, dim, ndvr, ndvr, runit, tnum1(1,1,m), ld1, pwfn, 1, czero, hpwfn, 1)
     call ZGBTRS('N', dim, ndvr, ndvr, 1, tden1(1,1,m), ld2, tpiv1(1,m), hpwfn, dim, info1)

     call ZGBMV('N', dim, dim, ndvr, ndvr, runit, tnum2(1,1,m), ld1, mwfn, 1, czero, hmwfn, 1)
     call ZGBTRS('N', dim, ndvr, ndvr, 1, tden2(1,1,m), ld2, tpiv2(1,m), hmwfn, dim, info2)

     wfn(1:dim, lval,     ifun) = half * (hpwfn(1:dim) + hmwfn(1:dim))
     wfn(1:dim, lval + 1, ifun) = half * (hpwfn(1:dim) - hmwfn(1:dim))
  end do

  deallocate(tden2)
  deallocate(tnum2)
  deallocate(tpiv2)
  deallocate(tden1)
  deallocate(tnum1)
  deallocate(tpiv1)
  deallocate(hmwfn)
  deallocate(hpwfn)
  deallocate(mwfn)
  deallocate(pwfn)

end subroutine adi_laser_vgauge_mixp
!######################################################################
