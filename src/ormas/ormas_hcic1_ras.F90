!######################################################################
subroutine ormas_hcic1_ras(int1e, cic, hcic)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nact,ndetx,nelact

  implicit none
  complex(c_double_complex), intent(in) :: int1e(1:nact, 1:nact)
  complex(c_double_complex), intent(in) :: cic(1:ndetx)
  complex(c_double_complex), intent(inout) :: hcic(1:ndetx)

  if (nact == 0) return
  if (nelact(2) > 0) call ormas_hcic1_ras_bbp(int1e, cic, hcic)
  if (nelact(1) > 0) call ormas_hcic1_ras_aap(int1e, cic, hcic)

!DEBUG
!  call ormas_cic_print(hcic, "test.h1cic")
!DEBUG

end subroutine ormas_hcic1_ras
!######################################################################
subroutine ormas_hcic1_ras_bbp(int1e, cic, hcic)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : thrcic
  use mod_bas, only : mtot
  use mod_ormas, only : nact,det_allowed,ndetx,ntot_alph_beta
  use mod_ormas, only : ndist_alph,nstr_alph
  use mod_ormas, only : ndist_beta,nstr_beta,dist_str_beta
  use mod_ormas, only : nstr_dist_m_alph,llstr_dist_m_alph,mval_beta,mmin_alph,mmax_alph
  use mod_ormas, only : n1x_m_beta,map1x_m_beta,p1x_beta,h1x_beta,eq1x_beta,sgn1x_beta

  implicit none
  complex(c_double_complex), intent(in) :: int1e(1:nact, 1:nact)
  complex(c_double_complex), intent(in) :: cic(1:ndetx)
  complex(c_double_complex), intent(inout) :: hcic(1:ndetx)
  integer(c_int) :: istr,jstr,kstr,i1x_m,i1x,ifun,jfun,idist,jdist,kdist,lla,ula,mvala
  complex(c_double_complex) :: tmp

  !$omp parallel default(shared) private(mvala,idist,i1x,ifun,jfun,kstr,kdist,tmp,lla,ula)
  !$omp do
  do istr = 1, nstr_beta
     mvala = mtot-mval_beta(istr)
     if (mvala>mmax_alph .or. mvala<mmin_alph) cycle

     idist = dist_str_beta(1,istr)
     do i1x_m = 1, n1x_m_beta(0,istr)
        i1x = map1x_m_beta(i1x_m,0,istr)
        ifun = h1x_beta(i1x,istr)
        jfun = p1x_beta(i1x,istr)
        kstr = eq1x_beta(i1x,istr)
        kdist = dist_str_beta(1,kstr)
        tmp = sgn1x_beta(i1x,istr) * int1e(ifun,jfun)
        if (abs(tmp) < thrcic) cycle
        do jdist = 1, ndist_alph
           if (det_allowed(jdist,kdist) == 0 .or. &
               det_allowed(jdist,idist) == 0) cycle
           lla = llstr_dist_m_alph(jdist,mvala)
           ula =  nstr_dist_m_alph(jdist,mvala) + lla - 1
           do jstr = lla, ula
!              hcic(mapf_detx(jstr,istr)) = &
!              hcic(mapf_detx(jstr,istr)) + tmp * &
!               cic(mapf_detx(jstr,kstr))
              hcic(ntot_alph_beta(istr)+jstr) = &
              hcic(ntot_alph_beta(istr)+jstr) + tmp * &
               cic(ntot_alph_beta(kstr)+jstr)
           end do
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ormas_hcic1_ras_bbp
!######################################################################
subroutine ormas_hcic1_ras_aap(int1e, cic, hcic)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : thrcic
  use mod_bas, only : mtot,smul
  use mod_ormas, only : nact,det_allowed,nelact,ndetx,mapr_detx,ntot_alph_beta
  use mod_ormas, only : ndist_beta,nstr_beta
  use mod_ormas, only : ndist_alph,nstr_alph,dist_str_alph
  use mod_ormas, only : nstr_dist_m_beta,llstr_dist_m_beta,mval_alph,mmin_beta,mmax_beta
  use mod_ormas, only : n1x_m_alph,map1x_m_alph,p1x_alph,h1x_alph,eq1x_alph,sgn1x_alph

  implicit none
  complex(c_double_complex), intent(in) :: int1e(1:nact,1:nact)
  complex(c_double_complex), intent(in) :: cic(1:ndetx)
  complex(c_double_complex), intent(inout) :: hcic(1:ndetx)
  integer(c_int) :: istr,jstr,kstr,i1x_m,i1x,ifun,jfun,idist,jdist,kdist,llb,ulb,idet,mvalb
  complex(c_double_complex) :: tmp
  complex(c_double_complex), allocatable :: dcic(:)

  !##### restricted case #####
!  if (smul == 1 .and. nelact(1) == nelact(2)) then
  if (.false.) then
     allocate(dcic(1:ndetx))
     dcic = hcic
     !$omp parallel default(shared)
     !$omp do
     do idet = 1, ndetx
!        hcic(idet) = hcic(idet) + dcic(mapf_detx(mapr_detx(2,idet),mapr_detx(1,idet)))
        hcic(idet) = hcic(idet) + dcic(ntot_alph_beta(mapr_detx(1,idet))+mapr_detx(2,idet))
     end do
!old     do istr = 1, nstr_beta
!old        do jstr = 1, istr
!old           hcic(mapf_detx(jstr,istr)) = &
!old           hcic(mapf_detx(jstr,istr)) + &
!old           hcic(mapf_detx(istr,jstr))
!old        end do
!old     end do
     !$omp end do
!!     !$omp do
!!     do istr = 1, nstr_beta
!!        do jstr = istr + 1, nstr_alph
!!           hcic(mapf_detx(jstr,istr)) = &
!!           hcic(mapf_detx(istr,jstr))
!!        end do
!!     end do
!!     !$omp end do
     !$omp end parallel
     deallocate(dcic)
     return
  end if

  !$omp parallel default(shared) private(idist,i1x,ifun,jfun,kstr,kdist,tmp,llb,ulb,mvalb)
  !$omp do
  do istr = 1, nstr_alph
     mvalb = mtot-mval_alph(istr)
     if (mvalb>mmax_beta .or. mvalb<mmin_beta) cycle

     idist = dist_str_alph(1,istr)
     do i1x_m = 1, n1x_m_alph(0,istr)
        i1x = map1x_m_alph(i1x_m,0,istr)
        ifun = h1x_alph(i1x,istr)
        jfun = p1x_alph(i1x,istr)
        kstr = eq1x_alph(i1x,istr)
        kdist = dist_str_alph(1,kstr)
        tmp = sgn1x_alph(i1x,istr) * int1e(ifun,jfun)
        if (abs(tmp) < thrcic) cycle
        do jdist = 1, ndist_beta
           if (det_allowed(kdist,jdist) == 0 .or. &
               det_allowed(idist,jdist) == 0) cycle
           llb = llstr_dist_m_beta(jdist,mvalb)
           ulb =  nstr_dist_m_beta(jdist,mvalb) + llb - 1
           do jstr = llb, ulb
!              hcic(mapf_detx(istr,jstr)) = &
!              hcic(mapf_detx(istr,jstr)) + tmp * &
!               cic(mapf_detx(kstr,jstr))
              hcic(ntot_alph_beta(jstr)+istr) = &
              hcic(ntot_alph_beta(jstr)+istr) + tmp * &
               cic(ntot_alph_beta(jstr)+kstr)
           end do
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ormas_hcic1_ras_aap
!######################################################################
