!######################################################################
subroutine hprod_ppop1(lfield,den1,orb,zorb,pzorb,gorb,vel,acc)

  use,intrinsic :: iso_c_binding
  use mod_bas,only : mval
  use mod_sph,only : lmax1
  use mod_rad,only : nrad,nradfc,ecs_flag,irad_ecs
  use mod_const,only : zero,czero,ctwo,runit
  use mod_ormas,only : nact,nfcore,ncore,nfun
  use mod_hprod,only : torb,v2orb
!  torb            ! scratch
! v2orb            ! scratch

  implicit none
  real(c_double),intent(in) :: lfield(1:3,1:3)
  complex(c_double_complex),intent(in) :: den1(1:nact,1:nact)
  complex(c_double_complex),intent(in) :: orb(1:(nrad-1),0:lmax1,1:nfun)
  complex(c_double_complex),intent(in) :: zorb(1:(nrad-1),0:lmax1,1:nfun)
  complex(c_double_complex),intent(in) :: pzorb(1:(nrad-1),0:lmax1,1:nfun)
  complex(c_double_complex),intent(in) :: gorb(1:(nrad-1),0:lmax1,1:nfun)
  real(c_double),intent(inout) :: vel,acc

  complex(c_double_complex) :: vtmp,vval
  complex(c_double_complex) :: atmp,aval
  complex(c_double_complex) :: atmp2,aval2
  complex(c_double_complex),parameter :: cfac = -runit
  integer(c_int) :: llrf,ulrf,llrd,ulrd,ifun,jfun,iact,jact,l,irad
  integer(c_int) :: rdim

  torb = 0D0
  v2orb = 0D0
  call hprod_projpp(runit,lfield,orb,torb)
  !call hprod_projpp(runit,lfield,zorb,v2orb)
  !call hprod_zprod_all(cfac,torb,v2orb)
  call hprod_zprod_all(runit,torb,v2orb)
  call hprod_projpp(cfac,lfield,zorb,v2orb)

! Orimo_ECS
  if (ecs_flag == 1) then
     rdim = irad_ecs - 1
  else 
     rdim = nrad - 1
  end if
! Orimo_ECS

  vval = czero
  aval = czero
  aval2 = czero
  !$omp parallel default(shared) private(ifun,jfun,llrf,ulrf,llrd,ulrd,vtmp,atmp,atmp2) reduction(+:vval,aval,aval2)
  !###########################
  call util_omp_disp(1,nradfc,llrf,ulrf)
! Orimo_ECS
  call util_omp_disp(1,rdim,llrd,ulrd)
! Orimo_ECS

  vtmp = czero
  atmp = czero
  atmp2 = czero
  do ifun = 1,nfcore
     do l = abs(mval(ifun)),lmax1
        do irad = llrf,ulrf
           vtmp = vtmp + conjg(zorb(irad,l,ifun)) * torb(irad,l,ifun)
           atmp = atmp + conjg(pzorb(irad,l,ifun)) * torb(irad,l,ifun)
           atmp2 = atmp2 + conjg(gorb(irad,l,ifun)) * v2orb(irad,l,ifun)
        end do
     end do
  end do

  do ifun = nfcore + 1,ncore
     do l = abs(mval(ifun)),lmax1
        do irad = llrd,ulrd
           vtmp = vtmp + conjg(zorb(irad,l,ifun)) * torb(irad,l,ifun)
           atmp = atmp + conjg(pzorb(irad,l,ifun)) * torb(irad,l,ifun)
           atmp2 = atmp2 + conjg(gorb(irad,l,ifun)) * v2orb(irad,l,ifun)
        end do
     end do
  end do

  vtmp = vtmp * 2D0
  atmp = atmp * 2D0
  atmp2 = atmp2 * 2D0

  do iact = 1,nact
     ifun = ncore + iact
     do jact = 1,nact
        jfun = ncore + jact
        if (mval(ifun) == mval(jfun)) then
           do l = abs(mval(ifun)),lmax1
              do irad = llrd,ulrd
                 vtmp = vtmp + conjg(zorb(irad,l,ifun)) * torb(irad,l,jfun) * den1(jact,iact)
                 atmp = atmp + conjg(pzorb(irad,l,ifun)) * torb(irad,l,jfun) * den1(jact,iact)
                 atmp2 = atmp2 + conjg(gorb(irad,l,ifun)) * v2orb(irad,l,jfun) * den1(jact,iact)
              end do
           end do
        end if
     end do
  end do
  vval = vval + vtmp
  aval = aval + atmp
  aval2 = aval2 + atmp2
  !###########################
  !$omp end parallel

  vel = vel + aimag(vval)*2.0
  acc = acc + aimag(aval)*2.0 - dble(aval2)*2.0
  !BUG: acc = acc + aimag(aval)*2.0 + dble(aval2)*2.0

end subroutine hprod_ppop1
!######################################################################
