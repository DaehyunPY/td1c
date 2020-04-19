!######################################################################
subroutine hprod_pzprod_all(zfac, wfn, hwfn)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : nbas
  use mod_rad, only : nrad
  use mod_ormas, only : nfun, nfcore

  implicit none
  complex(c_double_complex), intent(in) :: zfac
  complex(c_double_complex), intent(in) :: wfn(1:nbas, 1:*)
  complex(c_double_complex), intent(inout) :: hwfn(1:nbas, 1:*)

  if (nfcore > 0) call hprod_pzprod_fc(zfac, wfn, hwfn)
  call hprod_pzprod_dyn(zfac, wfn, hwfn)

end subroutine hprod_pzprod_all
!######################################################################
subroutine hprod_pzprod_dyn(zfac, wfn, hwfn)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : nbas
  use mod_rad, only : nrad
  use mod_ormas, only : nfcore, nfun

  implicit none
  complex(c_double_complex), intent(in) :: zfac
  complex(c_double_complex), intent(in) :: wfn(1:nbas, 1:*)
  complex(c_double_complex), intent(inout) :: hwfn(1:nbas, 1:*)

  call hprod_pzprod(zfac, wfn, hwfn, nfcore + 1, nfun, 1, nrad - 1)

end subroutine hprod_pzprod_dyn
!######################################################################
subroutine hprod_pzprod_fc(zfac, wfn, hwfn)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : nbas
  use mod_rad, only : nradfc
  use mod_ormas, only : nfcore

  implicit none
  complex(c_double_complex), intent(in) :: zfac
  complex(c_double_complex), intent(in) :: wfn(1:nbas, 1:*)
  complex(c_double_complex), intent(inout) :: hwfn(1:nbas, 1:*)

  if (nfcore == 0) return
  call hprod_pzprod(zfac, wfn, hwfn, 1, nfcore, 1, nradfc)

end subroutine hprod_pzprod_fc
!######################################################################
subroutine hprod_pzprod_fc1(zfac, wfn, hwfn)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : nbas
  use mod_rad, only : nradfc
  use mod_ormas, only : nfcore1, nfcore2, nfcore

  implicit none
  complex(c_double_complex), intent(in) :: zfac
  complex(c_double_complex), intent(in) :: wfn(1:nbas, 1:*)
  complex(c_double_complex), intent(inout) :: hwfn(1:nbas, 1:*)

  if (nfcore1 == 0) return
  call hprod_pzprod(zfac, wfn, hwfn, nfcore2+1, nfcore, 1, nradfc)

end subroutine hprod_pzprod_fc1
!######################################################################
subroutine hprod_pzprod_fc2(zfac, wfn, hwfn)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : nbas
  use mod_rad, only : nradfc
  use mod_ormas, only : nfcore1, nfcore2, nfcore

  implicit none
  complex(c_double_complex), intent(in) :: zfac
  complex(c_double_complex), intent(in) :: wfn(1:nbas, 1:*)
  complex(c_double_complex), intent(inout) :: hwfn(1:nbas, 1:*)

  if (nfcore2 == 0) return
  call hprod_pzprod(zfac, wfn, hwfn, 1, nfcore2, 1, nradfc)

end subroutine hprod_pzprod_fc2
!######################################################################
subroutine hprod_pzprod(zfac, wfn, hwfn, llfun, ulfun, llr, ulr)

  use, intrinsic :: iso_c_binding
  use mod_sph, only : lmax1
  use mod_hprod, only : h1orb_out
  use mod_rad, only : nrad, ecs_flag
  use mod_control, only : ioorot

  implicit none
  integer(c_long), intent(in) :: llfun, ulfun, llr, ulr
  complex(c_double_complex), intent(in) :: zfac
  complex(c_double_complex), intent(in) :: wfn(1:(nrad-1), 0:lmax1, 1:*)
  complex(c_double_complex), intent(inout) :: hwfn(1:(nrad-1), 0:lmax1, 1:*)
  integer(c_long) :: llrp, ulrp

  !$omp parallel default(shared) private(llrp, ulrp)
  !###########################
  call util_omp_disp(llr, ulr, llrp, ulrp)
  call hprod_pzprodp(zfac, wfn, hwfn, llfun, ulfun, llrp, ulrp)
  !###########################
  !$omp end parallel

  if (ecs_flag==1 .and. ioorot.ne.1) then
     call hprod_pzprod_ecsout(zfac, wfn, h1orb_out, llfun, ulfun)
  end if

end subroutine hprod_pzprod
!######################################################################
subroutine hprod_pzprodp(zfac, wfn, hwfn, llfun, ulfun, llrp, ulrp)

  use, intrinsic :: iso_c_binding
  use mod_sph, only : lmax1
  use mod_const, only : czero
  use mod_bas, only : mval, pmat, bas_pzfac1, bas_pzfac2
  use mod_rad, only : nrad, ndvr, xrad, wrad, cxrad, ecs_flag

  implicit none
  integer(c_long), intent(in) :: llfun, ulfun, llrp, ulrp
  complex(c_double_complex), intent(in) :: zfac
  complex(c_double_complex), intent(in) :: wfn(1:(nrad-1), 0:lmax1, 1:*)
  complex(c_double_complex), intent(inout) :: hwfn(1:(nrad-1), 0:lmax1, 1:*)

  integer(c_long) :: ifun, l, m, irad, jrad, jll, jul, jb1, kb1
  complex(c_double_complex) :: tmp1, tmp2, dwc
  complex(c_double_complex), allocatable :: dwfn(:,:)
  complex(c_double_complex), external :: zdotu

  allocate(dwfn(llrp:ulrp, 0:lmax1))

  if(ecs_flag == 0) then
     do ifun = llfun, ulfun
        m = mval(ifun)
        ! radial first derivative
        dwfn(llrp:ulrp, 0:lmax1) = czero
        do l = abs(m), lmax1
           do irad = llrp, ulrp
              kb1 = ndvr + 1 - irad
              jll = max(1,        irad - ndvr)
              jul = min(nrad - 1, irad + ndvr)
!new
              dwfn(irad, l) = dwfn(irad, l) + zdotu(jul-jll+1, wfn(jll,l,ifun), 1, pmat(jll+kb1,irad), 1)
!def              dwc = czero
!def              do jrad = jll, jul
!def                 jb1 = kb1 + jrad
!def                 dwc = dwc + wfn(jrad, l, ifun) * pmat(jb1, irad)
!def              end do
!def              dwfn(irad, l) = dwfn(irad, l) + dwc
!new
           end do
        end do
        ! summation
        do l = abs(m), lmax1 - 1
           tmp1 = zfac * bas_pzfac1(l, m)
           do irad = llrp, ulrp
              tmp2 = zfac * bas_pzfac2(irad, l, m)
              hwfn(irad, l,     ifun) = hwfn(irad, l,     ifun) &
                                      + dwfn(irad, l + 1)       * tmp1 &
                                      +  wfn(irad, l + 1, ifun) * tmp2
              hwfn(irad, l + 1, ifun) = hwfn(irad, l + 1, ifun) &
                                      + dwfn(irad, l)           * tmp1 &
                                      -  wfn(irad, l,     ifun) * tmp2
           end do
        end do
     end do
! ORIMO_ECS
  else if(ecs_flag == 1) then
  !else if(0) then
     ! There is no difference from the above.
     do ifun = llfun, ulfun
        m = mval(ifun)
        ! radial first derivative
        dwfn(llrp:ulrp, 0:lmax1) = czero
        do l = abs(m), lmax1
           do irad = llrp, ulrp
              kb1 = ndvr + 1 - irad
              jll = max(1,        irad - ndvr)
              jul = min(nrad - 1, irad + ndvr)
!new
              dwfn(irad, l) = dwfn(irad, l) + zdotu(jul-jll+1, wfn(jll,l,ifun), 1, pmat(jll+kb1,irad), 1)
!def              dwc = czero
!def              do jrad = jll, jul
!def                 jb1 = kb1 + jrad
!def                 dwc = dwc + wfn(jrad, l, ifun) * pmat(jb1, irad)
!def              end do
!def              dwfn(irad, l) = dwfn(irad, l) + dwc
!new
           end do
        end do
        ! summation
        do l = abs(m), lmax1 - 1
           tmp1 = zfac * bas_pzfac1(l, m)
           do irad = llrp, ulrp
              tmp2 = zfac * bas_pzfac2(irad, l, m)
              hwfn(irad, l,     ifun) = hwfn(irad, l,     ifun) &
                                      + dwfn(irad, l + 1)       * tmp1 &
                                      +  wfn(irad, l + 1, ifun) * tmp2
              hwfn(irad, l + 1, ifun) = hwfn(irad, l + 1, ifun) &
                                      + dwfn(irad, l)           * tmp1 &
                                      -  wfn(irad, l,     ifun) * tmp2
           end do
        end do
     end do
! ORIMO_ECS
  end if

  deallocate(dwfn)
  
end subroutine hprod_pzprodp
!######################################################################
subroutine hprod_pzprod_ecsout(zfac, wfn, hwfn, llfun, ulfun)

  use, intrinsic :: iso_c_binding
  use mod_sph, only : lmax1
  use mod_bas, only : mval, pmat, bas_pzfac1
  use mod_rad, only : nrad, ndvr, xrad, wrad, cxrad, ecs_flag, irad_ecs

  implicit none
  integer(c_long), intent(in) :: llfun, ulfun
  complex(c_double_complex), intent(in) :: zfac
  complex(c_double_complex), intent(in) :: wfn(1:(nrad-1), 0:lmax1, 1:*)
  complex(c_double_complex), intent(out) :: hwfn(1:(nrad-1), 0:lmax1, 1:*)

  integer(c_long) :: ifun, l, irad, jrad, lll, ull
  complex(c_double_complex) :: tmp1
  complex(c_double_complex), allocatable :: dwfn(:,:,:)

  allocate(dwfn(irad_ecs-ndvr:irad_ecs-1,0:lmax1,llfun:ulfun))
  dwfn = 0d0

  !$omp parallel default(shared) private(tmp1)
  !$omp do
  !###########################
  do ifun = llfun, ulfun
     do l = abs(mval(ifun)), lmax1
        do irad = irad_ecs-ndvr, irad_ecs-1
           do jrad = irad_ecs, irad+ndvr
              dwfn(irad,l,ifun) = dwfn(irad,l,ifun) + wfn(jrad,l,ifun)*pmat(ndvr+1+jrad-irad,irad)
           end do
        end do
     end do
     hwfn(:,:,ifun) = 0d0
     do l = abs(mval(ifun)), lmax1-1
        tmp1 = zfac*bas_pzfac1(l,mval(ifun))
        do irad = irad_ecs-ndvr, irad_ecs-1
           hwfn(irad,l,  ifun) = hwfn(irad,l,  ifun) + dwfn(irad,l+1,ifun)*tmp1
           hwfn(irad,l+1,ifun) = hwfn(irad,l+1,ifun) + dwfn(irad,l,  ifun)*tmp1
        end do
     end do
  end do
  !###########################
  !$omp end do
  !$omp end parallel

  deallocate(dwfn)
  
end subroutine hprod_pzprod_ecsout
!######################################################################
