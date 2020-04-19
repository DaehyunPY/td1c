!######################################################################
subroutine adi_laser_lgauge(icomp, istag, dtime, lfield, wfn)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad
  use mod_sph, only : lmax1
  use mod_ormas, only : nfun

  implicit none
  integer(c_int), intent(in) :: icomp, istag
  real(c_double), intent(in) :: dtime, lfield(1:3)
  complex(c_double_complex), intent(inout) :: wfn(1:(nrad-1), 0:lmax1, 1:nfun)

  integer(c_int) :: eve_odd, ll_eo, ul_eo, pp_eo
  integer(c_int) :: lmaxh, lllo2, ullo2, lo2, lval

  lmaxh = lmax1 / 2
  if (istag == 0) then
     ll_eo = 0
     ul_eo = 1
     pp_eo =+1
  else
     ll_eo = 1
     ul_eo = 0
     pp_eo =-1
  end if

!DEBUG
!  write(6, "('adi_laser_lgauge:', 2i10)") lmax1, lmaxh
!  do eve_odd = 1, 0, -1
!     do lo2 = 0, lmaxh
!        lval = eve_odd + 2 * lo2
!        if (lval >= lmax1) cycle
!        write(6, "('     ', 3i10)") eve_odd, lo2, lval
!     end do
!  end do  
!  stop "Stop for debug @ adi_laser_lgauge."
!DEBUG

  do eve_odd = ll_eo, ul_eo, pp_eo
     !$omp parallel default(shared) private(lval, lllo2, ullo2)
     !###########################
     call util_omp_disp(0, lmaxh, lllo2, ullo2)
     do lo2 = lllo2, ullo2
        lval = eve_odd + 2 * lo2
        if (lval >= lmax1) cycle
        call adi_laser_lgaugep(icomp, dtime, lfield, wfn, lval)
     end do
     !###########################
     !$omp end parallel
  end do

end subroutine adi_laser_lgauge
!######################################################################
subroutine adi_laser_lgaugep(icomp, dtime, lfield, wfn, lval)

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

  facd = runit / sqrt((two * lval + one) * (two * lval + three))
  if (icomp == 1) then
     facl = iunit * lfield(3) * dtime * half * facd
  else
     facl = runit * lfield(3) * dtime * half * facd
  end if

  do ifun = nfcore + 1, nfun
     if (lval < abs(mval(ifun))) cycle
     beta = facl * sqrt((lval + one) ** two - mval(ifun) ** two)
     do irad = 1, nrad - 1
        tmp = beta * xrad(irad)
        ombd = one - tmp
        opbd = one + tmp
        gammp = ombd / opbd
        gammm = opbd / ombd

        wfn0 = wfn(irad, lval,     ifun)
        wfn1 = wfn(irad, lval + 1, ifun)
        yp = (wfn0 + wfn1) * gammp
        ym = (wfn0 - wfn1) * gammm
        wfn(irad, lval,     ifun) = half * (yp + ym)
        wfn(irad, lval + 1, ifun) = half * (yp - ym)
     end do
  end do

end subroutine adi_laser_lgaugep
!######################################################################
