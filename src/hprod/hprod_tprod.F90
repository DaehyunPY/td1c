!######################################################################
subroutine hprod_tprod_all(wfn, hwfn)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : nbas
  use mod_rad, only : nrad
  use mod_ormas, only : nfun, nfcore

  implicit none
  complex(c_double_complex), intent(in) :: wfn(1:nbas, 1:*)
  complex(c_double_complex), intent(inout) :: hwfn(1:nbas, 1:*)

! call hprod_tprod(wfn, hwfn, 1, nfun, 1, nrad - 1)
 if (nfcore > 0) call hprod_tprod_fc(wfn, hwfn)
 call hprod_tprod_dyn(wfn, hwfn)

end subroutine hprod_tprod_all
!######################################################################
subroutine hprod_tprod_dyn(wfn, hwfn)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : nbas
  use mod_rad, only : nrad
  use mod_ormas, only : nfcore, nfun

  implicit none
  complex(c_double_complex), intent(in) :: wfn(1:nbas, 1:*)
  complex(c_double_complex), intent(inout) :: hwfn(1:nbas, 1:*)

  call hprod_tprod(wfn, hwfn, nfcore + 1, nfun, 1, nrad - 1)

end subroutine hprod_tprod_dyn
!######################################################################
subroutine hprod_tprod_fc(wfn, hwfn)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : nbas
  use mod_rad, only : nradfc
  use mod_ormas, only : nfcore

  implicit none
  complex(c_double_complex), intent(in) :: wfn(1:nbas, 1:*)
  complex(c_double_complex), intent(inout) :: hwfn(1:nbas, 1:*)

  if (nfcore == 0) return
  call hprod_tprod(wfn, hwfn, 1, nfcore, 1, nradfc)

end subroutine hprod_tprod_fc
!######################################################################
subroutine hprod_tprod_fc1(wfn, hwfn)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : nbas
  use mod_rad, only : nradfc
  use mod_ormas, only : nfcore1, nfcore2, nfcore

  implicit none
  complex(c_double_complex), intent(in) :: wfn(1:nbas, 1:*)
  complex(c_double_complex), intent(inout) :: hwfn(1:nbas, 1:*)

  if (nfcore1 == 0) return
  call hprod_tprod(wfn, hwfn, nfcore2+1, nfcore, 1, nradfc)

end subroutine hprod_tprod_fc1
!######################################################################
subroutine hprod_tprod_fc2(wfn, hwfn)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : nbas
  use mod_rad, only : nradfc
  use mod_ormas, only : nfcore1, nfcore2, nfcore

  implicit none
  complex(c_double_complex), intent(in) :: wfn(1:nbas, 1:*)
  complex(c_double_complex), intent(inout) :: hwfn(1:nbas, 1:*)

  if (nfcore2 == 0) return
  call hprod_tprod(wfn, hwfn, 1, nfcore2, 1, nradfc)

end subroutine hprod_tprod_fc2
!######################################################################
subroutine hprod_tprod_ene(wfn, hwfn)

  use, intrinsic :: iso_c_binding
  use mod_sph, only : lmax1
  use mod_hprod, only : v2jfc
  use mod_rad, only : nrad, nradfc
  use mod_ormas, only : nfun, nfcore
  use mod_control, only : jfc_implicit

  implicit none
  complex(c_double_complex), intent(in) :: wfn(1:(nrad-1), 0:lmax1, 1:*)
  complex(c_double_complex), intent(inout) :: hwfn(1:(nrad-1), 0:lmax1, 1:*)
  integer(c_long) :: llrf, ulrf, llrd, ulrd, ifun, l, irad

! call hprod_tprod(wfn, hwfn, 1, nfun, 1, nrad - 1)
 if (nfcore > 0) call hprod_tprod_fc(wfn, hwfn)
 call hprod_tprod_dyn(wfn, hwfn)

 if (.not. jfc_implicit) return

 !$omp parallel default(shared) private(llrf, ulrf, llrd, ulrd)
 !###########################
 call util_omp_disp(1, nradfc,   llrf, ulrf)
 call util_omp_disp(1, nrad - 1, llrd, ulrd)
 do ifun = 1, nfcore
    do l = 0, lmax1
       do irad = llrf, ulrf
          hwfn(irad, l, ifun) = hwfn(irad, l, ifun) - v2jfc(irad) * wfn(irad, l, ifun)
       end do
    end do
 end do
 do ifun = nfcore + 1, nfun
    do l = 0, lmax1
       do irad = llrd, ulrd
          hwfn(irad, l, ifun) = hwfn(irad, l, ifun) - v2jfc(irad) * wfn(irad, l, ifun)
       end do
    end do
 end do
 !###########################
 !$omp end parallel

end subroutine hprod_tprod_ene
!######################################################################
subroutine hprod_tprod(wfn, hwfn, llfun, ulfun, llr, ulr)

  use, intrinsic :: iso_c_binding
  use mod_sph, only : lmax1
  use mod_hprod, only : h0orb_out
  use mod_control, only : ioorot
  use mod_rad, only : nrad, ecs_flag

  implicit none
  integer(c_long), intent(in) :: llfun, ulfun, llr, ulr
  complex(c_double_complex), intent(in) :: wfn(1:(nrad-1), 0:lmax1, 1:*)
  complex(c_double_complex), intent(inout) :: hwfn(1:(nrad-1), 0:lmax1, 1:*)
  integer(c_long) :: llrp, ulrp

  !$omp parallel default(shared) private(llrp, ulrp)
  !###########################
  call util_omp_disp(llr, ulr, llrp, ulrp)
  call hprod_tprodp(wfn, hwfn, llfun, ulfun, llrp, ulrp)
  !###########################
  !$omp end parallel

  if (ecs_flag==1 .and. ioorot==0) then
     call hprod_tprod_ecsout(wfn, h0orb_out, llfun, ulfun)
  end if

end subroutine hprod_tprod
!######################################################################
subroutine hprod_tprodp(wfn, hwfn, llfun, ulfun, llrp, ulrp)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad, ndvr
  use mod_bas, only : mval, tmat
  use mod_sph, only : lmax1
  use mod_const, only : czero
  use mod_rad, only : mapf, mapb

  implicit none
  integer(c_long), intent(in) :: llfun, ulfun, llrp, ulrp
  complex(c_double_complex), intent(in) :: wfn(1:(nrad-1), 0:lmax1, 1:*)
  complex(c_double_complex), intent(inout) :: hwfn(1:(nrad-1), 0:lmax1, 1:*)

  complex(c_double_complex) :: hwc
  integer(c_long) :: ifun, irad, jrad, jb1, kb1, l, jll, jul
  integer(c_long) :: jb, ife, iloc, jloc
  complex(c_double_complex), external :: zdotu

!check here  do ifun = llfun, ulfun
!check here     do l = abs(mval(ifun)), lmax1
!check here        do irad = llrp, ulrp
!check here           kb1 = ndvr + 1 - irad
!check here           jll = max(1,        irad - ndvr)
!check here           jul = min(nrad - 1, irad + ndvr)
!check here
!check here           hwc = czero
!check here           do jrad = jll, jul
!check here              jb1 = kb1 + jrad
!check here              hwc = hwc + wfn(jrad, l, ifun) * tmat(jb1, irad, l)
!check here           end do
!check here           hwfn(irad, l, ifun) = hwfn(irad, l, ifun) + hwc
!check here        end do
!check here     end do
!check here  end do

!new
!Why hwfn is INITIALIZED?
  do ifun = llfun, ulfun
     do l = abs(mval(ifun)), lmax1
        do irad = llrp, ulrp
           kb1 = ndvr + 1 - irad
           jll = max(1,        irad - ndvr)
           jul = min(nrad - 1, irad + ndvr)
           hwfn(irad, l, ifun) = zdotu(jul-jll+1, wfn(jll,l,ifun), 1, tmat(jll+kb1,irad,l), 1)
        end do
     end do
  end do
!def  do ifun = llfun, ulfun
!def     do l = abs(mval(ifun)), lmax1
!def        do irad = llrp, ulrp
!def           ife = mapb(irad)
!def           iloc = irad - mapf(ife);
!def
!def           hwc = czero
!def
!def           ! non bridge elements
!def           do jloc = 0, ndvr
!def              jrad = mapf(ife) + jloc;
!def              if (jrad > 0 .and. jrad < nrad) then !NEWNEWNEW
!def                 jb = ndvr + 1 + jrad - irad
!def                 hwc = hwc + wfn(jrad, l, ifun) * tmat(jb, irad, l)
!def              end if
!def           end do
!def
!def           ! bridge functions
!def           if (iloc == 0 .and. ife /= 0) then
!def              !do jloc = 0, ndvr
!def              do jloc = 0, ndvr - 1
!def                 jrad = mapf(ife - 1) + jloc;
!def                 if (jrad > 0 .and. jrad < nrad) then !NEWNEWNEW
!def                    jb = ndvr + 1 + jrad - irad
!def                    hwc = hwc + wfn(jrad, l, ifun) * tmat(jb, irad, l)
!def                 end if
!def              end do
!def           end if
!def
!def           hwfn(irad, l, ifun) = hwc
!def        end do
!def     end do
!def  end do
!new

end subroutine hprod_tprodp
!######################################################################
subroutine hprod_tprod_ecsout(wfn, hwfn, llfun, ulfun)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad, ndvr, irad_ecs, ndvr
  use mod_bas, only : mval, tmat
  use mod_sph, only : lmax1

  implicit none
  integer(c_long), intent(in) :: llfun, ulfun
  complex(c_double_complex), intent(in) :: wfn(1:(nrad-1), 0:lmax1, 1:*)
  complex(c_double_complex), intent(out) :: hwfn(1:(nrad-1), 0:lmax1, 1:*)
  integer(c_long) :: ifun, irad, jrad, l, lll, ull

  !$omp parallel default(shared) private(lll,ull)
  !###########################
  call util_omp_disp(0, lmax1, lll, ull)
  do ifun = llfun, ulfun
     do l = lll, ull
        hwfn(:,l,ifun) = 0d0
        if (l < mval(ifun)) cycle
        do irad = irad_ecs-ndvr, irad_ecs-1
           do jrad = irad_ecs, irad+ndvr
              hwfn(irad,l,ifun) = hwfn(irad,l,ifun) + wfn(jrad,l,ifun)*tmat(ndvr+1+jrad-irad,irad,l)
           end do
        end do
     end do
  end do
  !###########################
  !$omp end parallel

end subroutine hprod_tprod_ecsout
!######################################################################
