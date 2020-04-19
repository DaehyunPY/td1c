!######################################################################
subroutine hprod_op1tr_init(nfun_,orb,dP_,vP_,aP_,dQ_,vQ_,aQ_)

  use,intrinsic :: iso_c_binding
  use mod_ormas,only : nfun,neltot
  use mod_sph,only : nlat,lmax1
  use mod_op1tr

  implicit none
  integer(c_int),target,intent(in) :: nfun_
  complex(c_double_complex),intent(in) :: orb(*)
  complex(c_double_complex),target,intent(in) :: dP_(1:nfun_,1:nfun_)
  complex(c_double_complex),target,intent(in) :: vP_(1:nfun_,1:nfun_)
  complex(c_double_complex),target,intent(in) :: aP_(1:nfun_,1:nfun_)
  complex(c_double_complex),target,intent(in) :: dQ_(1:nfun_)
  complex(c_double_complex),target,intent(in) :: vQ_(1:nfun_)
  complex(c_double_complex),target,intent(in) :: aQ_(1:nfun_)

  do_op1tr = .true.
  op1tr_nfun => nfun_
  op1tr_dP => dP_
  op1tr_vP => vP_
  op1tr_aP => aP_
  op1tr_dQ => dQ_
  op1tr_vQ => vQ_
  op1tr_aQ => aQ_

  call hprod_op1tr_nrad(orb)

  allocate(op1tr_orb0(1:op1tr_nrad,0:lmax1,1:nfun))
  allocate(op1tr_orb (1:op1tr_nrad,0:lmax1,1:nfun))
  allocate(op1tr_dorb(1:op1tr_nrad,0:lmax1,1:nfun))
  allocate(op1tr_vorb(1:op1tr_nrad,0:lmax1,1:nfun))
  allocate(op1tr_aorb(1:op1tr_nrad,0:lmax1,1:nfun))
  allocate(op1tr_orbg (1:op1tr_nrad,1:nlat,1:nfun))
  allocate(op1tr_dorbg(1:op1tr_nrad,1:nlat,1:nfun))
  allocate(op1tr_vorbg(1:op1tr_nrad,1:nlat,1:nfun))
  allocate(op1tr_aorbg(1:op1tr_nrad,1:nlat,1:nfun))
  allocate(op1tr_umat(1:op1tr_nfun,1:nfun))
  allocate(op1tr_dmat(1:op1tr_nfun,1:nfun))
  allocate(op1tr_vmat(1:op1tr_nfun,1:nfun))
  allocate(op1tr_amat(1:op1tr_nfun,1:nfun))
  allocate(op1tr_dmat0(1:op1tr_nfun,1:op1tr_nfun))
  allocate(op1tr_vmat0(1:op1tr_nfun,1:op1tr_nfun))
  allocate(op1tr_amat0(1:op1tr_nfun,1:op1tr_nfun))
  allocate(den1_tr0(1:nfun,1:nfun))
  allocate(den1_tr1(1:nfun,1:op1tr_nfun))
  allocate(den1_tr2(1:op1tr_nfun,1:op1tr_nfun))

  call hprod_op1tr_opmat0()

end subroutine hprod_op1tr_init
!######################################################################
subroutine hprod_op1tr_final()

  use,intrinsic :: iso_c_binding
  use mod_op1tr

  implicit none

  deallocate(op1tr_orb0)
  deallocate(op1tr_orb )
  deallocate(op1tr_dorb)
  deallocate(op1tr_vorb)
  deallocate(op1tr_aorb)
  deallocate(op1tr_orbg )
  deallocate(op1tr_dorbg)
  deallocate(op1tr_vorbg)
  deallocate(op1tr_aorbg)
  deallocate(op1tr_umat)
  deallocate(op1tr_dmat)
  deallocate(op1tr_vmat)
  deallocate(op1tr_amat)
  deallocate(op1tr_dmat0)
  deallocate(op1tr_vmat0)
  deallocate(op1tr_amat0)
  deallocate(den1_tr0)
  deallocate(den1_tr1)
  deallocate(den1_tr2)
!  deallocate(op1tr_dP)
!  deallocate(op1tr_vP)
!  deallocate(op1tr_aP)
!  deallocate(op1tr_dQ)
!  deallocate(op1tr_vQ)
!  deallocate(op1tr_aQ)

end subroutine hprod_op1tr_final
!######################################################################
subroutine hprod_op1tr(lfield,dip,vel,acc)

  use,intrinsic :: iso_c_binding
  use mod_const,only : iunit
  use mod_rad,only : xrad,nrad
  use mod_sph,only : lmax1,nlat,cost,legf1,legb1
  use mod_bas,only : mval,lval
  use mod_control,only : igauge
  use mod_ormas,only : nfcore,ncore,nfun,nact,neltot
  use mod_op1tr

  implicit none
  real(c_double), intent(in) :: lfield(1:3, 1:3)
  real(c_double),intent(inout) :: dip(1:3),vel(1:3),acc(1:3)

  integer(c_int) :: ifun,jfun,pfun,qfun,mi,li,iact,jact
  integer(c_int) :: llr,ulr,irad,ilat,l,lmax_init
  complex(c_double_complex) :: faclv,corb0

  op1tr_umat(:,:) = 0d0
  op1tr_dmat(:,:) = 0d0
  op1tr_vmat(:,:) = 0d0
  op1tr_amat(:,:) = 0d0
  do ifun = 1,op1tr_nfun
     mi = mval(ifun)
     li = lval(ifun)
     do pfun = 1,nfun
        if (mval(pfun) .ne. mi) cycle
        do irad = 1,op1tr_nrad
           corb0 = conjg(op1tr_orb0(irad,li,ifun))
           op1tr_umat(ifun,pfun) = op1tr_umat(ifun,pfun) + corb0 * op1tr_orb (irad,li,pfun)
           op1tr_dmat(ifun,pfun) = op1tr_dmat(ifun,pfun) + corb0 * op1tr_dorb(irad,li,pfun)
           op1tr_vmat(ifun,pfun) = op1tr_vmat(ifun,pfun) + corb0 * op1tr_vorb(irad,li,pfun)
           op1tr_amat(ifun,pfun) = op1tr_amat(ifun,pfun) + corb0 * op1tr_aorb(irad,li,pfun)
        end do
     end do
  end do

!debug
!  do ifun = 1,op1tr_nfun
!     write(6,"('# umat:  ',i5)",advance='no') ifun
!     do pfun = 1,nfun
!        write(6,"(f15.5)",advance='no') dble(op1tr_umat(ifun,pfun))
!     end do
!     write(6,*)
!  end do
!  write(6,*)
!
!  do ifun = 1,op1tr_nfun
!     write(6,"('# amat:  ',i5)",advance='no') ifun
!     do pfun = 1,nfun
!        write(6,"(f15.5)",advance='no') dble(op1tr_amat(ifun,pfun))
!     end do
!     write(6,*)
!  end do
!  stop
!debug

  den1_tr1(:,:) = 0d0
  do ifun = 1,op1tr_nfun
     mi = mval(ifun)
     do pfun = 1,nfun
        if (mval(pfun) .ne. mi) cycle
        do qfun = 1,nfun
           if (mval(qfun) .ne. mi) cycle
           den1_tr1(pfun,ifun) = den1_tr1(pfun,ifun) + den1_tr0(pfun,qfun)*conjg(op1tr_umat(ifun,qfun))
        end do
     end do
  end do

  den1_tr2(:,:) = 0d0
  do ifun = 1,op1tr_nfun
     mi = mval(ifun)
     do jfun = 1,op1tr_nfun
        if (mval(jfun) .ne. mi) cycle
        do pfun = 1,nfun
           if (mval(pfun) .ne. mi) cycle
           den1_tr2(jfun,ifun) = den1_tr2(jfun,ifun) + op1tr_umat(jfun,pfun)*den1_tr1(pfun,ifun)
        end do
     end do
  end do

  op1tr_dP(:,:) = 0d0
  op1tr_vP(:,:) = 0d0
  op1tr_aP(:,:) = 0d0
  do ifun = 1,op1tr_nfun
     li = lval(ifun)
     mi = mval(ifun)
     do jfun = 1,op1tr_nfun
        if (lval(jfun).eq.li .or. mval(jfun).ne.mi) cycle
        op1tr_dP(ifun,jfun) = op1tr_dmat0(ifun,jfun)*den1_tr2(jfun,ifun)
        op1tr_vP(ifun,jfun) = op1tr_vmat0(ifun,jfun)*den1_tr2(jfun,ifun)
        op1tr_aP(ifun,jfun) = op1tr_amat0(ifun,jfun)*den1_tr2(jfun,ifun)
     end do
     !op1tr_aP(ifun,ifun) = op1tr_aP(ifun,ifun) - lfield(3,2)*den1_tr2(ifun,ifun)
  end do

  op1tr_dQ(:) = 0d0
  op1tr_vQ(:) = 0d0
  op1tr_aQ(:) = 0d0  
  do ifun = 1,op1tr_nfun
     li = lval(ifun)
     mi = mval(ifun)
     do pfun = 1,nfun
        if (mval(pfun).ne.mi) cycle
        op1tr_dQ(ifun) = op1tr_dQ(ifun) + op1tr_dmat(ifun,pfun)*den1_tr1(pfun,ifun)
        op1tr_vQ(ifun) = op1tr_vQ(ifun) + op1tr_vmat(ifun,pfun)*den1_tr1(pfun,ifun)
        op1tr_aQ(ifun) = op1tr_aQ(ifun) + op1tr_amat(ifun,pfun)*den1_tr1(pfun,ifun)
     end do
     do jfun = 1,op1tr_nfun
        if (lval(jfun).eq.li .or. mval(jfun).ne.mi) cycle
        op1tr_dQ(ifun) = op1tr_dQ(ifun) - op1tr_dP(ifun,jfun)
        op1tr_vQ(ifun) = op1tr_vQ(ifun) - op1tr_vP(ifun,jfun)
        op1tr_aQ(ifun) = op1tr_aQ(ifun) - op1tr_aP(ifun,jfun)
     end do
  end do

!  do ifun = 1,op1tr_nfun
!     op1tr_aP(ifun,ifun) = op1tr_aP(ifun,ifun) - lfield(3,2)*den1_tr2(ifun,ifun)
!  end do

  !summation
!  dip(1) = 0d0
!  vel(1) = 0d0
!  acc(1) = -lfield(3, 2)*neltot(3)
  dip(1) = 0d0
  vel(1) = 0d0
  acc(1) = 0d0
  do ifun = 1,op1tr_nfun
     !li = lval(ifun)
     mi = mval(ifun)
     do jfun = 1,op1tr_nfun
        !if (lval(jfun).eq.li .or. mval(jfun).ne.mi) cycle
        if (mval(jfun).ne.mi) cycle
        dip(1) = dip(1) + dble(op1tr_dP(ifun,jfun))
        vel(1) = vel(1) + dble(op1tr_vP(ifun,jfun))
        acc(1) = acc(1) + dble(op1tr_aP(ifun,jfun))
     end do
     dip(1) = dip(1) + 2d0*dble(op1tr_dQ(ifun))
     vel(1) = vel(1) + 2d0*dble(op1tr_vQ(ifun))
     acc(1) = acc(1) + 2d0*dble(op1tr_aQ(ifun))
  end do

  dip(3) = dip(1) + dip(2)
  vel(3) = vel(1) + vel(2)
  acc(3) = acc(1) + acc(2)

end subroutine hprod_op1tr
!######################################################################
subroutine hprod_op1tr_orbin(lfield,den1,orb,dorb,vorb,aorb)

  use,intrinsic :: iso_c_binding
  use mod_const,only : iunit
  use mod_rad,only : xrad,nrad
  use mod_sph,only : lmax1,nlat,cost,legf1,legb1
  use mod_bas,only : mval,lval
  use mod_control,only : igauge
  use mod_ormas,only : nfcore,ncore,nfun,nact
  use mod_op1tr

  implicit none
  real(c_double),intent(in) :: lfield(1:3,1:3)
  complex(c_double_complex),intent(in) :: den1(1:nact,1:nact)
  complex(c_double_complex),intent(in) :: orb(1:(nrad-1),0:lmax1,1:nfun)
  complex(c_double_complex),intent(in) :: dorb(1:(nrad-1),0:lmax1,1:nfun)
  complex(c_double_complex),intent(in) :: vorb(1:(nrad-1),0:lmax1,1:nfun)
  complex(c_double_complex),intent(in) :: aorb(1:(nrad-1),0:lmax1,1:nfun)

  integer(c_int) :: ifun,jfun,pfun,qfun,mi,li,iact,jact
  integer(c_int) :: llr,ulr,irad,ilat,l,lmax_init
  complex(c_double_complex) :: faclv,corb0

  den1_tr0(:,:) = 0d0
  do ifun = 1,ncore
     den1_tr0(ifun,ifun) = 2d0
  end do
  do iact = 1,nact
     ifun = ncore + iact
     mi = mval(ifun)
     do jact = 1,nact
        jfun = ncore + jact
        if (mval(jfun) .ne. mi) cycle
        den1_tr0(jfun,ifun) = den1(jact,iact)
     end do
  end do

  op1tr_orb  = 0d0
  op1tr_dorb = 0d0
  op1tr_vorb = 0d0
  op1tr_aorb = 0d0
  do ifun = 1,nfun
     mi = mval(ifun)
     do l = abs(mi),lmax1
        do irad = 1,op1tr_nrad
           op1tr_orb (irad,l,ifun) = orb (irad,l,ifun)
           op1tr_dorb(irad,l,ifun) = dorb(irad,l,ifun)
           op1tr_vorb(irad,l,ifun) = vorb(irad,l,ifun)
           op1tr_aorb(irad,l,ifun) = aorb(irad,l,ifun)
        end do
     end do
  end do

  if (igauge == 1) then
     op1tr_orbg  = 0d0
     op1tr_dorbg = 0d0
     op1tr_vorbg = 0d0
     op1tr_aorbg = 0d0
     !$omp parallel default(shared) private(mi,llr,ulr)
     call util_omp_disp(1,op1tr_nrad,llr,ulr)
     do ifun = 1,nfun
        mi = mval(ifun)
        do ilat = 1,nlat
           do l = abs(mi),lmax1
              do irad = llr,ulr
                 op1tr_orbg (irad,ilat,ifun) = op1tr_orbg (irad,ilat,ifun) + op1tr_orb (irad,l,ifun) * legf1(l,ilat,mi)
                 op1tr_dorbg(irad,ilat,ifun) = op1tr_dorbg(irad,ilat,ifun) + op1tr_dorb(irad,l,ifun) * legf1(l,ilat,mi)
                 op1tr_vorbg(irad,ilat,ifun) = op1tr_vorbg(irad,ilat,ifun) + op1tr_vorb(irad,l,ifun) * legf1(l,ilat,mi)
                 op1tr_aorbg(irad,ilat,ifun) = op1tr_aorbg(irad,ilat,ifun) + op1tr_aorb(irad,l,ifun) * legf1(l,ilat,mi)
              end do
           end do
        end do
     end do
     !$omp end parallel
     
     !$omp parallel default(shared) private(faclv,llr,ulr)
     !###########################
     call util_omp_disp(1,op1tr_nrad,llr,ulr)
     do ilat = 1,nlat
        do irad = llr,ulr
           faclv = exp(iunit * lfield(3,3) * xrad(irad) * cost(ilat))
           do ifun = 1,nfun
              op1tr_orbg (irad,ilat,ifun) = op1tr_orbg (irad,ilat,ifun) * faclv
              op1tr_dorbg(irad,ilat,ifun) = op1tr_dorbg(irad,ilat,ifun) * faclv
              op1tr_vorbg(irad,ilat,ifun) = op1tr_vorbg(irad,ilat,ifun) * faclv
              op1tr_aorbg(irad,ilat,ifun) = op1tr_aorbg(irad,ilat,ifun) * faclv
           end do
        end do
     end do
     !###########################
     !$omp end parallel
     
     op1tr_orb  = 0d0
     op1tr_dorb = 0d0
     op1tr_vorb = 0d0
     op1tr_aorb = 0d0
     !$omp parallel default(shared) private(lmax_init,mi,llr,ulr)
     lmax_init = maxval(lval(1:nfun))
     call util_omp_disp(1,op1tr_nrad,llr,ulr)
     do ifun = 1,nfun
        mi = mval(ifun)
        do l = abs(mi),lmax_init
           do ilat = 1,nlat
              do irad = llr,ulr
                 op1tr_orb (irad,l,ifun) = op1tr_orb (irad,l,ifun) + op1tr_orbg (irad,ilat,ifun) * legb1(ilat,l,mi)
                 op1tr_dorb(irad,l,ifun) = op1tr_dorb(irad,l,ifun) + op1tr_dorbg(irad,ilat,ifun) * legb1(ilat,l,mi)
                 op1tr_vorb(irad,l,ifun) = op1tr_vorb(irad,l,ifun) + op1tr_vorbg(irad,ilat,ifun) * legb1(ilat,l,mi)
                 op1tr_aorb(irad,l,ifun) = op1tr_aorb(irad,l,ifun) + op1tr_aorbg(irad,ilat,ifun) * legb1(ilat,l,mi)
              end do
           end do
        end do
     end do
     !$omp end parallel
  end if

end subroutine hprod_op1tr_orbin
!######################################################################
subroutine hprod_op1tr_opmat0()

  use,intrinsic :: iso_c_binding
  use mod_const,only : runit
  use mod_bas,only : mval,lval,nval
  use mod_hprod,only : orb0,h0orb,h1orb,v2orb
  use mod_op1tr

  implicit none
  integer(c_int) :: ifun,jfun,irad,li,lj
  complex(c_double_complex) :: corb0

  op1tr_orb0 = orb0
  h0orb = 0d0
  h1orb = 0d0
  v2orb = 0d0
  call hprod_zprod_all(runit,orb0,h0orb)
  call hprod_pzprod_all(runit,orb0,h1orb)
  call hprod_azprod_all(runit,orb0,v2orb)

  op1tr_dmat0 = 0d0
  op1tr_vmat0 = 0d0
  op1tr_amat0 = 0d0
  do ifun = 1,op1tr_nfun
     li = lval(ifun)
     do jfun = 1,op1tr_nfun
        lj = lval(jfun)
        if (lj.eq.li .or. mval(jfun).ne.mval(ifun)) cycle
        do irad = 1,op1tr_nrad
           corb0 = conjg(op1tr_orb0(irad,li,ifun))
           op1tr_dmat0(ifun,jfun) = op1tr_dmat0(ifun,jfun) + corb0 * h0orb(irad,li,jfun)
           op1tr_vmat0(ifun,jfun) = op1tr_vmat0(ifun,jfun) + corb0 * h1orb(irad,li,jfun)
           op1tr_amat0(ifun,jfun) = op1tr_amat0(ifun,jfun) + corb0 * v2orb(irad,li,jfun)
        end do
     end do
  end do

!debug
!  do ifun = 1,op1tr_nfun
!     write(6,"('# amat0: ',i5)",advance='no') ifun
!     do jfun = 1,op1tr_nfun
!        write(6,"(f15.5)",advance='no') dble(op1tr_amat0(ifun,jfun))
!     end do
!     write(6,*)
!  end do
!debug

!debug
!  write(6,"('hprod_op1tr_opmat0: dmat0')")
!  do ifun = 1,op1tr_nfun
!     do jfun = 1,op1tr_nfun
!        write(6,"(4i5,' :',4i5,2f20.10)") &
!             ifun,nval(ifun),lval(ifun),mval(ifun), &
!             jfun,nval(jfun),lval(jfun),mval(jfun), &
!             op1tr_dmat0(ifun,jfun)
!     end do
!  end do
!  write(6,"('hprod_op1tr_opmat0: vmat0')")
!  do ifun = 1,op1tr_nfun
!     do jfun = 1,op1tr_nfun
!        write(6,"(4i5,' :',4i5,2f20.10)") &
!             ifun,nval(ifun),lval(ifun),mval(ifun), &
!             jfun,nval(jfun),lval(jfun),mval(jfun), &
!             op1tr_vmat0(ifun,jfun)
!     end do
!  end do
!  write(6,"('hprod_op1tr_opmat0: amat0')")
!  do ifun = 1,op1tr_nfun
!     do jfun = 1,op1tr_nfun
!        write(6,"(4i5,' :',4i5,2f20.10)") &
!             ifun,nval(ifun),lval(ifun),mval(ifun), &
!             jfun,nval(jfun),lval(jfun),mval(jfun), &
!             op1tr_amat0(ifun,jfun)
!     end do
!  end do
!debug

end subroutine hprod_op1tr_opmat0
!######################################################################
subroutine hprod_op1tr_nrad(orb)

  use,intrinsic :: iso_c_binding
  use mod_rad,only : nrad,xrad,mapf,mapb,ndvr,ecs_flag,irad_ecs
  use mod_const,only : zero
  use mod_sph,only : lmax1
  use mod_bas,only : mval
  use mod_op1tr,only : op1tr_nrad,op1tr_nfun

  implicit none
  complex(c_double_complex),intent(in) :: orb(1:(nrad-1),0:lmax1,1:*)

  real(c_double) :: vtail
  real(c_double),parameter :: tiny = 1.0D-8
!  real(c_double),parameter :: tiny = 1.0D-18
! real(c_double),parameter :: thrtail = 1.0D-20
! real(c_double),parameter :: thrtail = 1.0D-15
! real(c_double),parameter :: thrtail = 1.0D-12
! real(c_double),parameter :: thrtail = 1.0D-10

  integer(c_int) :: iorb,irad,l,nradmax,lradmax,nradtmp,max_irad

! Sato_ECS
  if (ecs_flag == 0) then
     max_irad = nrad - 1
  else
     max_irad = irad_ecs - 1
  end if
  op1tr_nrad = max_irad
! Sato_ECS

  do iorb = 1,op1tr_nfun
     nradmax = 1
     lradmax = abs(mval(iorb))
     do l = abs(mval(iorb)),lmax1
        nradtmp = max_irad
        do irad = 1,max_irad
           vtail = sqrt(conjg(orb(irad,l,iorb)) * orb(irad,l,iorb))
           !debug
           !write(6,"('hprod_op1tr_nrad: ',3i5,f20.10,2e20.10)") iorb,l,irad,xrad(irad),vtail,tiny
           !debug
           if (vtail < tiny) then
              nradtmp = irad
              exit
           end if
        end do
        if (nradtmp > nradmax) then
           nradmax = nradtmp
           lradmax = l
        end if
     end do
     op1tr_nrad = min(max_irad,nradmax)
  end do
  op1tr_nrad = min(max_irad,mapf(mapb(op1tr_nrad)) + ndvr)
  write(6,"('hprod_op1tr_nrad: op1tr_nrad = ',i5,F20.10,E20.10,E20.10)") op1tr_nrad,xrad(op1tr_nrad),vtail,tiny
  
end subroutine hprod_op1tr_nrad
!######################################################################
