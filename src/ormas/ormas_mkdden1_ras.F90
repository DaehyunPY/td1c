!######################################################################
subroutine ormas_mkdden1_ras(fac, cic, dcic, den1)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : smul
  use mod_const, only : czero
  use mod_ormas, only : ndetx
  use mod_ormas, only : ncore,nsub,nelact,nact,lorb_sub,mval

  implicit none
  complex(c_double_complex), intent(in) :: fac
  complex(c_double_complex), intent(in) :: cic(1:ndetx)
  complex(c_double_complex), intent(in) :: dcic(1:ndetx)
  complex(c_double_complex), intent(inout) :: den1(1:nact, 1:nact)

  integer(c_long) :: iproc,nproc,iact,jact,isub,jsub,istr,jstr
  integer(c_long), external :: util_omp_nproc
  complex(c_double_complex) , allocatable :: den1p(:,:,:)

  if (nact == 0) return
  nproc = util_omp_nproc()
  allocate(den1p(1:nact, 1:nact, 0:(nproc-1)))

  den1p = czero
  if (nelact(2) > 0) call ormas_mkdden1_ras_bb(cic, dcic, den1p, nproc)
  if (smul==1 .and. nelact(1)==nelact(2)) then
!  if (.false.) then
     den1p = den1p + den1p
  else
     if (nelact(1) > 0) call ormas_mkdden1_ras_aa(cic, dcic, den1p, nproc)
  end if

  do iproc = 1, nproc - 1
     den1p(:,:,0) = den1p(:,:,0) + den1p(:,:,iproc)
  end do

  ! unify the bra-derivative contributions, which are a.h.c. of ket-derivatives
  do isub = 1, nsub
     do jsub = 1, isub - 1
        do iact = lorb_sub(1,isub), lorb_sub(2,isub)
           do jact = lorb_sub(1,jsub), lorb_sub(2,jsub)
              den1p(jact,iact,0) = den1p(jact,iact,0) - conjg(den1p(iact,jact,0))
              den1p(iact,jact,0) = -conjg(den1p(jact,iact,0))
           end do
        end do
        do iact = lorb_sub(1,isub), lorb_sub(2,isub)
           do jact = lorb_sub(1,jsub), lorb_sub(2,jsub)
              if (mval(ncore+iact) == mval(ncore+jact)) then
                 den1(iact,jact) = den1(iact,jact) + den1p(iact,jact,0) * fac
                 den1(jact,iact) = den1(jact,iact) + den1p(jact,iact,0) * fac
              end if
           end do
        end do
     end do
  end do

  deallocate(den1p)

end subroutine ormas_mkdden1_ras
!######################################################################
!######################################################################
subroutine ormas_mkdden1_ras_bb(cic, dcic, den1p, nproc)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : mtot
  use mod_ormas, only : ndetx,ntot_alph_beta
  use mod_ormas, only : nact,det_allowed,ndist_alph
  use mod_ormas, only : nstr_dist_m_alph,llstr_dist_m_alph,dist_str_beta
  use mod_ormas, only : nstr_alph,nstr_beta,p1x_beta,h1x_beta,eq1x_beta,sgn1x_beta
  use mod_ormas, only : n1x_m_beta,map1x_m_beta,mval_beta,mmin_alph,mmax_alph

  implicit none
  integer(c_long), intent(in) :: nproc
  complex(c_double_complex) , intent(in) :: cic(1:ndetx)
  complex(c_double_complex) , intent(in) :: dcic(1:ndetx)
  complex(c_double_complex) , intent(inout) :: den1p(1:nact,1:nact,0:(nproc-1))

  integer(c_long) :: istr,jstr,kstr,i1x,i1x_m,iact,jact,iproc
  integer(c_long) :: idist,jdist,kdist,mvala,lla,ula
  integer(c_long), external :: util_omp_iproc

  !$omp parallel default(shared) private(mvala,idist,i1x,iact,jact,kstr,kdist,lla,ula,iproc)
  iproc = util_omp_iproc()
  !$omp do 
  do istr = 1, nstr_beta
     mvala = mtot-mval_beta(istr)
     if (mvala>mmax_alph .or. mvala<mmin_alph) cycle

     idist = dist_str_beta(1,istr)
     do i1x_m = 1, n1x_m_beta(0,istr)
        i1x = map1x_m_beta(i1x_m,0,istr)
        iact = h1x_beta (i1x,istr)
        jact = p1x_beta (i1x,istr)
        kstr = eq1x_beta(i1x,istr)
        kdist = dist_str_beta(1,kstr)
        do jdist = 1, ndist_alph
           if (det_allowed(jdist,idist) == 0 .or. &
               det_allowed(jdist,kdist) == 0) cycle
           lla = llstr_dist_m_alph(jdist,mvala)
           ula =  nstr_dist_m_alph(jdist,mvala) + lla - 1
           do jstr = lla, ula
!              den1p(jact,iact,iproc) = &
!              den1p(jact,iact,iproc) + conjg(cic(mapf_detx(jstr,istr))) &
!                   * sgn1x_beta(i1x,istr) * dcic(mapf_detx(jstr,kstr))
              den1p(jact,iact,iproc) = &
              den1p(jact,iact,iproc) + conjg(cic(ntot_alph_beta(istr)+jstr)) &
                   * sgn1x_beta(i1x,istr) * dcic(ntot_alph_beta(kstr)+jstr)
           end do
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ormas_mkdden1_ras_bb
!######################################################################
subroutine ormas_mkdden1_ras_aa(cic, dcic, den1p, nproc)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : mtot
  use mod_ormas, only : ndetx,ntot_alph_beta
  use mod_ormas, only : nact,det_allowed,ndist_beta
  use mod_ormas, only : nstr_dist_m_beta,llstr_dist_m_beta,dist_str_alph
  use mod_ormas, only : nstr_beta,nstr_alph,p1x_alph,h1x_alph,eq1x_alph,sgn1x_alph
  use mod_ormas, only : n1x_m_alph,map1x_m_alph,mval_alph,mmin_beta,mmax_beta

  implicit none
  integer(c_long), intent(in) :: nproc
  complex(c_double_complex) , intent(in) :: cic(1:ndetx)
  complex(c_double_complex) , intent(in) :: dcic(1:ndetx)
  complex(c_double_complex) , intent(inout) :: den1p(1:nact,1:nact,0:(nproc-1))

  integer(c_long) :: istr,jstr,kstr,i1x,i1x_m,iact,jact,iproc
  integer(c_long) :: idist,jdist,kdist,dista,mvalb,llb,ulb
  integer(c_long), external :: util_omp_iproc

  !$omp parallel default(shared) private(mvalb,idist,i1x,iact,jact,kstr,kdist,llb,ulb,iproc)
  iproc = util_omp_iproc()
  !$omp do 
  do istr = 1, nstr_alph
     mvalb = mtot-mval_alph(istr)
     if (mvalb>mmax_beta .or. mvalb<mmin_beta) cycle

     idist = dist_str_alph(1,istr)
     do i1x_m = 1, n1x_m_alph(0,istr)
        i1x = map1x_m_alph(i1x_m,0,istr)
        iact = h1x_alph (i1x,istr)
        jact = p1x_alph (i1x,istr)
        kstr = eq1x_alph(i1x,istr)
        kdist = dist_str_alph(1,kstr)
        do jdist = 1, ndist_beta
           if (det_allowed(idist,jdist) == 0 .or. &
               det_allowed(kdist,jdist) == 0) cycle
           llb = llstr_dist_m_beta(jdist,mvalb)
           ulb =  nstr_dist_m_beta(jdist,mvalb) + llb - 1
           do jstr = llb, ulb
!              den1p(jact,iact,iproc) = &
!              den1p(jact,iact,iproc) + conjg(cic(mapf_detx(istr,jstr))) &
!                   * sgn1x_alph(i1x,istr) * dcic(mapf_detx(kstr,jstr))
              den1p(jact,iact,iproc) = &
              den1p(jact,iact,iproc) + conjg(cic(ntot_alph_beta(jstr)+istr)) &
                   * sgn1x_alph(i1x,istr) * dcic(ntot_alph_beta(jstr)+kstr)
           end do
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ormas_mkdden1_ras_aa
!######################################################################
