!######################################################################
subroutine ormas_mkden1_ras(cic, den1)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : smul
  use mod_const, only : czero
  use mod_ormas, only : nelact,nact

  implicit none
  complex(c_double_complex) , intent(in) :: cic(1:*)
  complex(c_double_complex) , intent(out) :: den1(1:nact, 1:nact)

  integer(c_int) :: iproc, nproc, iact, jact
  integer(c_int), external :: util_omp_nproc
  complex(c_double_complex) , allocatable :: den1p(:,:,:)

  if (nact == 0) return
  nproc = util_omp_nproc()
  allocate(den1p(1:nact, 1:nact, 0:(nproc-1)))
  call util_zcopy(nproc*nact**2, czero, 0, den1p, 1)

  if (nelact(2) >= 1) call ormas_mkden1_ras_bb(cic, den1p, nproc)

  if (smul==1 .and. nelact(1)==nelact(2)) then
     den1p = den1p + den1p
  else
     if (nelact(1) >= 1) call ormas_mkden1_ras_aa(cic, den1p, nproc)
  end if

  den1 = czero
  do iproc = 0, nproc - 1
     den1(:,:) = den1(:,:) + den1p(:,:,iproc)
  end do
  deallocate(den1p)

!DEBUG
!  call util_print_vec(nact**2, den1, "test.den1")
!DEBUG

end subroutine ormas_mkden1_ras
!######################################################################
!######################################################################
subroutine ormas_mkden1_ras_bb(cic, den1p, nproc)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : mtot
  use mod_ormas, only : nact,det_allowed,ndist_alph,ndetx,ntot_alph_beta
  use mod_ormas, only : nstr_dist_m_alph,llstr_dist_m_alph,dist_str_beta
  use mod_ormas, only : nstr_alph,nstr_beta,p1x_beta,h1x_beta,eq1x_beta,sgn1x_beta
  use mod_ormas, only : n1x_m_beta,map1x_m_beta,mval_beta,mmin_alph,mmax_alph

  implicit none
  integer(c_int), intent(in) :: nproc
  complex(c_double_complex) , intent(in) :: cic(1:ndetx)
  complex(c_double_complex) , intent(inout) :: den1p(1:nact,1:nact,0:(nproc-1))

  integer(c_int) :: istr,jstr,kstr,i1x,i1x_m,iact,jact,iproc
  integer(c_int) :: jdist,lla,ula,mvala
  integer(c_int), external :: util_omp_iproc

!write(6,"('nstr_beta:  ', i5)") nstr_beta
!do istr = 1, nstr_beta
!   write(6,"('n1x_m_beta: ', 2i5)") istr,n1x_m_beta(0,istr)
!end do
!do istr = 1, nstr_beta
!   do i1x_m = 1, n1x_m_beta(0,istr)
!      i1x = map1x_m_beta(i1x_m,0,istr)
!      iact = h1x_beta (i1x,istr)
!      jact = p1x_beta (i1x,istr)
!      kstr = eq1x_beta(i1x,istr)
!      write(6,"('i1x_m_beta: ', 6i5)") istr,i1x_m,i1x,iact,jact,kstr
!   end do
!end do
!do istr = 1, nstr_beta
!   mvala = mtot-mval_beta(istr)
!   if (mvala > mmax_alph .or. mvala < mmin_alph) cycle
!   do i1x_m = 1, n1x_m_beta(0,istr)
!      i1x = map1x_m_beta(i1x_m,0,istr)
!      iact = h1x_beta (i1x,istr)
!      jact = p1x_beta (i1x,istr)
!      kstr = eq1x_beta(i1x,istr)
!      do jdist = 1, ndist_alph
!         if (det_allowed(jdist,dist_str_beta(1,istr)) == 0 .or. &
!             det_allowed(jdist,dist_str_beta(1,kstr)) == 0) cycle
!         lla = llstr_dist_m_alph(jdist,mvala)
!         ula =  nstr_dist_m_alph(jdist,mvala) + lla - 1
!         write(6,"('range of istra for this istrb/kstrb: ', 8i5)") istr,kstr,jdist,mtot,mval_beta(istr),mvala, &
!              llstr_dist_m_alph(jdist,mvala), &
!              nstr_dist_m_alph(jdist,mvala)
!      end do
!   end do
!end do

  !$omp parallel default(shared) private(i1x,iact,jact,kstr,lla,ula,mvala,iproc)
  iproc = util_omp_iproc()
  !$omp do 
  do istr = 1, nstr_beta
     mvala = mtot-mval_beta(istr)
     if (mvala > mmax_alph .or. mvala < mmin_alph) cycle

     do i1x_m = 1, n1x_m_beta(0,istr)
        i1x = map1x_m_beta(i1x_m,0,istr)
        iact = h1x_beta (i1x,istr)
        jact = p1x_beta (i1x,istr)
        kstr = eq1x_beta(i1x,istr)
        do jdist = 1, ndist_alph
           if (det_allowed(jdist,dist_str_beta(1,istr)) == 0 .or. &
               det_allowed(jdist,dist_str_beta(1,kstr)) == 0) cycle
           lla = llstr_dist_m_alph(jdist,mvala)
           ula =  nstr_dist_m_alph(jdist,mvala) + lla - 1
           do jstr = lla, ula
!              den1p(jact,iact,iproc) = &
!              den1p(jact,iact,iproc) + conjg(cic(mapf_detx(jstr,istr))) &
!                    * sgn1x_beta(i1x,istr) * cic(mapf_detx(jstr,kstr))
              den1p(jact,iact,iproc) = &
              den1p(jact,iact,iproc) + conjg(cic(ntot_alph_beta(istr)+jstr)) &
                    * sgn1x_beta(i1x,istr) * cic(ntot_alph_beta(kstr)+jstr)
!write(6,"(2i5,':',3i5,f20.10,3i5,f20.10,i5)") &
!     iact,jact, &
!     ntot_alph_beta(istr)+jstr,jstr,istr,dble(cic(ntot_alph_beta(istr)+jstr)), &
!     ntot_alph_beta(kstr)+jstr,jstr,kstr,dble(cic(ntot_alph_beta(kstr)+jstr)), &
!     sgn1x_beta(i1x,istr)
           end do
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ormas_mkden1_ras_bb
!######################################################################
subroutine ormas_mkden1_ras_aa(cic, den1p, nproc)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : mtot
  use mod_ormas, only : nact,det_allowed,ndist_beta,ndetx,ntot_alph_beta
  use mod_ormas, only : nstr_dist_m_beta,llstr_dist_m_beta,dist_str_alph
  use mod_ormas, only : nstr_beta,nstr_alph,p1x_alph,h1x_alph,eq1x_alph,sgn1x_alph
  use mod_ormas, only : n1x_m_alph,map1x_m_alph,mval_alph,mmin_beta,mmax_beta

  implicit none
  integer(c_int), intent(in) :: nproc
  complex(c_double_complex) , intent(in) :: cic(1:ndetx)
  complex(c_double_complex) , intent(inout) :: den1p(1:nact,1:nact,0:(nproc-1))

  integer(c_int) :: istr,jstr,kstr,i1x,i1x_m,iact,jact,iproc
  integer(c_int) :: idist,jdist,kdist,mvalb,llb,ulb
  integer(c_int), external :: util_omp_iproc

  !$omp parallel default(shared) private(mvalb,idist,i1x,iact,jact,kstr,kdist,llb,ulb,iproc)
  iproc = util_omp_iproc()
  !$omp do 
  do istr = 1, nstr_alph
     mvalb = mtot-mval_alph(istr)
     if (mvalb > mmax_beta .or. mvalb < mmin_beta) cycle

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
!                    * sgn1x_alph(i1x,istr) * cic(mapf_detx(kstr,jstr))
              den1p(jact,iact,iproc) = &
              den1p(jact,iact,iproc) + conjg(cic(ntot_alph_beta(jstr)+istr)) &
                    * sgn1x_alph(i1x,istr) * cic(ntot_alph_beta(jstr)+kstr)
           end do
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ormas_mkden1_ras_aa
!######################################################################
