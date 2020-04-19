!#######################################################################
subroutine hprod_dipole_ipx(max_ipx,rad_ipx,lfield,wfn,cic,dip,vel,acc)

  use,intrinsic :: iso_c_binding
  use mod_bas,only : nbas,mval
  use mod_ormas,only : neltot,nfun,nfcore,ndcore
  use mod_const,only : zero,two,four,czero,runit
  use mod_control,only : igauge,fedvr_normalized
  use mod_hprod,only : den1,orb,orbg,h0orb,h1orb,gorb,torb,v2orb,ovlp

  implicit none
  integer(c_int),intent(in) :: max_ipx
  real(c_double),intent(in) :: rad_ipx
  real(c_double),intent(in) :: lfield(1:3,1:3)
  complex(c_double_complex),intent(in) :: wfn(*)
  complex(c_double_complex),intent(in) :: cic(*)
  complex(c_double_complex),external :: hprod_trace1
  real(c_double),intent(out) :: dip(0:max_ipx),vel(0:max_ipx),acc(0:max_ipx)

  integer(c_int) :: ifun,jfun,iipx
  real(c_double) :: ipx(0:max_ipx)
  complex(c_double) :: denipx(1:nfun,1:nfun,0:max_ipx)
  complex(c_double) :: dtmp,dmat(1:nfun,1:nfun)
  complex(c_double) :: vtmp,vmat(1:nfun,1:nfun)
  complex(c_double) :: atmp,amat(1:nfun,1:nfun)

  if (nfcore > 0) then
     write(6,"('WARNING: hprod_dipole_ipx with FCs')")
  end if

  call hprod_orbin(lfield,wfn,orb,orbg)
  call hprod_mkovlp(rad_ipx,orb,orb,ovlp);
  call ormas_denipx(max_ipx,ovlp,cic,ipx,denipx)

!debug
!  call ormas_mkden1(cic, den1)
!  do ifun = 1,nfun
!  do jfun = 1,nfun
!     write(6,"('hprod_dipole_ipx: ',2i5)",advance='no') ifun,jfun
!     write(6,"(2f12.8)",advance='no') ovlp(ifun,jfun)
!     write(6,"(2f12.8)",advance='no') den1(ifun,jfun)
!     write(6,"(2f12.8)",advance='no') &
!          denipx(ifun,jfun,0)+denipx(ifun,jfun,1)+ &
!          denipx(ifun,jfun,2)+denipx(ifun,jfun,3)+denipx(ifun,jfun,4)
!!     do iipx = 0, max_ipx
!!        write(6,"(2f12.8)",advance='no') denipx(ifun,jfun,iipx)
!!     end do
!     write(6,*)
!  end do
!  end do
!  stop
!debug

!old  if (nfcore > 0) then
!old     stop "hprod_dipole_ipx: nyi for nfcore > 0."
!old  else if (ndcore > 0) then
!old     stop "hprod_dipole_ipx: nyi for ndcore > 0."
!old  else
!old     call ormas_mkden1_ipx(ovlp,cic,den1);
!old  end if

  ! property vectors
  call zclear_omp(nbas*nfun,h0orb)
  call zclear_omp(nbas*nfun,h1orb)
  call zclear_omp(nbas*nfun,v2orb)
  call hprod_zprod_all(runit,orb,h0orb)
  call hprod_pzprod_all(runit,orb,h1orb)
  call hprod_azprod_all(runit,orb,v2orb)

  call hprod_mkovlp(rad_ipx,orb,h0orb,dmat);
  call hprod_mkovlp(rad_ipx,orb,h1orb,vmat);
  call hprod_mkovlp(rad_ipx,orb,v2orb,amat);

  do iipx = 0, max_ipx
     dtmp = 0.0
     vtmp = 0.0
     atmp = 0.0
     do ifun = 1, nfun
        do jfun = 1, nfun
           if (mval(ifun) .ne. mval(jfun)) cycle
           dtmp = dtmp + dmat(ifun,jfun)*denipx(jfun,ifun,iipx)
           vtmp = vtmp + vmat(ifun,jfun)*denipx(jfun,ifun,iipx)
           atmp = atmp + amat(ifun,jfun)*denipx(jfun,ifun,iipx)
        end do
     end do
     dip(iipx) = dble(dtmp)
     vel(iipx) = dble(vtmp)
     acc(iipx) = dble(atmp)
! The laser contribution is neglected since they exhibit
! spurious charge-resolved contributions to non-linear response
! which are summed up to zero.
!     if (igauge == 1) &
!     vel(iipx) = vel(iipx) + lfield(3,3)*neltot(3)*ipx(iipx)
!     acc(iipx) = acc(iipx) - lfield(3,2)*neltot(3)*ipx(iipx)
  end do

end subroutine hprod_dipole_ipx
!#######################################################################
