!#######################################################################
subroutine hprod_dipole_ipd(max_ipd,rad_ipd,lfield,wfn,cic,dip,vel,acc)

  use,intrinsic :: iso_c_binding
  use mod_bas,only : nbas,mval
  use mod_ormas,only : neltot,nfun,nfcore,ndcore
  use mod_const,only : zero,two,four,czero,runit
  use mod_control,only : igauge,fedvr_normalized
  use mod_hprod,only : den1,orb,orbg,h0orb,h1orb,gorb,torb,v2orb,ovlp

  implicit none
  integer(c_int),intent(in) :: max_ipd
  real(c_double),intent(in) :: rad_ipd
  real(c_double),intent(in) :: lfield(1:3,1:3)
  complex(c_double_complex),intent(in) :: wfn(*)
  complex(c_double_complex),intent(in) :: cic(*)
  complex(c_double_complex),external :: hprod_trace1
  real(c_double),intent(out) :: dip(0:max_ipd),vel(0:max_ipd),acc(0:max_ipd)

  integer(c_int) :: ifun,jfun,iipd
  real(c_double) :: ipd(0:max_ipd)
  complex(c_double) :: denipd(1:nfun,1:nfun,0:max_ipd)
  complex(c_double) :: dtmp,dmat(1:nfun,1:nfun)
  complex(c_double) :: vtmp,vmat(1:nfun,1:nfun)
  complex(c_double) :: atmp,amat(1:nfun,1:nfun)

  if (nfcore > 0) then
     write(6,"('WARNING: hprod_dipole_ipd with FCs')")
  end if

  call hprod_orbin(lfield,wfn,orb,orbg)
  call hprod_mkovlp(rad_ipd,orb,orb,ovlp);
  call ormas_denipd(max_ipd,ovlp,cic,ipd,denipd)

!debug
!  call ormas_mkden1(cic, den1)
!  do ifun = 1,nfun
!  do jfun = 1,nfun
!     write(6,"('hprod_dipole_ipd: ',2i5)",advance='no') ifun,jfun
!     write(6,"(2f12.8)",advance='no') ovlp(ifun,jfun)
!     write(6,"(2f12.8)",advance='no') den1(ifun,jfun)
!     write(6,"(2f12.8)",advance='no') &
!          denipd(ifun,jfun,0)+denipd(ifun,jfun,1)+ &
!          denipd(ifun,jfun,2)+denipd(ifun,jfun,3)+denipd(ifun,jfun,4)
!!     do iipd = 0, max_ipd
!!        write(6,"(2f12.8)",advance='no') denipd(ifun,jfun,iipd)
!!     end do
!     write(6,*)
!  end do
!  end do
!  stop
!debug

!old  if (nfcore > 0) then
!old     stop "hprod_dipole_ipd: nyi for nfcore > 0."
!old  else if (ndcore > 0) then
!old     stop "hprod_dipole_ipd: nyi for ndcore > 0."
!old  else
!old     call ormas_mkden1_ipd(ovlp,cic,den1);
!old  end if

  ! property vectors
  call zclear_omp(nbas*nfun,h0orb)
  call zclear_omp(nbas*nfun,h1orb)
  call zclear_omp(nbas*nfun,v2orb)
  call hprod_zprod_all(runit,orb,h0orb)
  call hprod_pzprod_all(runit,orb,h1orb)
  call hprod_azprod_all(runit,orb,v2orb)

  call hprod_mkovlp(rad_ipd,orb,h0orb,dmat);
  call hprod_mkovlp(rad_ipd,orb,h1orb,vmat);
  call hprod_mkovlp(rad_ipd,orb,v2orb,amat);


! The laser contribution is included...
  do ifun = 1, nfun
     if (igauge == 1) &
     vmat(ifun,ifun) = vmat(ifun,ifun) + lfield(3,3)
     amat(ifun,ifun) = amat(ifun,ifun) - lfield(3,2)
  end do

  do iipd = 0, max_ipd
     dtmp = 0.0
     vtmp = 0.0
     atmp = 0.0
     do ifun = 1, nfun
        do jfun = 1, nfun
           if (mval(ifun) .ne. mval(jfun)) cycle
           dtmp = dtmp + dmat(ifun,jfun)*denipd(jfun,ifun,iipd)
           vtmp = vtmp + vmat(ifun,jfun)*denipd(jfun,ifun,iipd)
           atmp = atmp + amat(ifun,jfun)*denipd(jfun,ifun,iipd)
        end do
     end do
     dip(iipd) = dble(dtmp)
     vel(iipd) = dble(vtmp)
     acc(iipd) = dble(atmp)
! The laser contribution is neglected since they exhibit
! spurious charge-resolved contributions to non-linear response
! which are summed up to zero.
!     if (igauge == 1) &
!     vel(iipd) = vel(iipd) + lfield(3,3)*neltot(3)*ipd(iipd)
!     acc(iipd) = acc(iipd) - lfield(3,2)*neltot(3)*ipd(iipd)
  end do

end subroutine hprod_dipole_ipd
!#######################################################################
