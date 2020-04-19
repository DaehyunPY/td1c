!######################################################################
subroutine ormas_mkden2x_ras(cic, den2)

  use, intrinsic :: iso_c_binding
  use mod_const, only : czero, ctwo
  use mod_ormas, only : s2zero, nact,nelact, nstr_alph, nstr_beta, ndetx

  implicit none
  complex(c_double_complex) , intent(in) :: cic(1:ndetx)
  complex(c_double_complex) , intent(out) :: den2(1:nact, 1:nact, 1:nact, 1:nact)

  integer(c_long) :: lwork
  integer(c_long) :: nproc,iproc,iact,jact,kact,lact
  complex(c_double_complex) , allocatable :: work(:,:)
  complex(c_double_complex) , allocatable :: den2aa(:,:,:,:,:)
  complex(c_double_complex) , allocatable :: den2bb(:,:,:,:,:)
  complex(c_double_complex) , allocatable :: den2ab(:,:,:,:,:)
  integer(c_long), external :: util_omp_nproc

  if (nact == 0) return
  nproc = util_omp_nproc()
  lwork = max(nstr_alph, nstr_beta)
  allocate(work(lwork, 0:(nproc-1)))
  allocate(den2aa(1:nact,1:nact,1:nact,1:nact,0:(nproc-1)))
  allocate(den2bb(1:nact,1:nact,1:nact,1:nact,0:(nproc-1)))
  allocate(den2ab(1:nact,1:nact,1:nact,1:nact,0:(nproc-1)))

  den2 = czero
  den2bb = czero
  den2aa = czero
  den2ab = czero

  if (nelact(2) > 0) call ormas_mkden2x_ras_bbp(cic,lwork,work,den2bb,nproc)
  if (nelact(1) > 0) call ormas_mkden2x_ras_aap(cic,lwork,work,den2aa,nproc)
  if (nelact(1) > 0 .and. nelact(2) > 0) call ormas_mkden2x_ras_abp(cic,lwork,work,den2ab,nproc)
  call ormas_mkden2x_ras_sum(den2aa,den2bb,den2ab,den2,nproc)

  deallocate(den2ab)
  deallocate(den2bb)
  deallocate(den2aa)
  deallocate(work)

end subroutine ormas_mkden2x_ras
!######################################################################
!######################################################################
subroutine ormas_mkden2x_ras_bbp(cic, lwork, work, den2p, nproc)

  use, intrinsic :: iso_c_binding
  use mod_sph, only : mmax2
  use mod_bas, only : mtot, smul
  use mod_ormas, only : ndetx,ntot_alph_beta
  use mod_ormas, only : nact,det_allowed,nstr_alph,nstr_beta,nelact
  use mod_ormas, only : ndist_alph,nstr_dist_m_alph,llstr_dist_m_alph
  use mod_ormas, only : dist_str_beta,n1x_m_beta,map1x_m_beta,mval_beta
  use mod_ormas, only : p1x_beta,h1x_beta,eq1x_beta,sgn1x_beta,mmin_alph,mmax_alph

  implicit none
  integer(c_long), intent(in) :: lwork, nproc
  complex(c_double_complex) , intent(in) :: cic(1:ndetx)
  complex(c_double_complex) , intent(out) :: work(1:lwork, 0:(nproc-1))
  complex(c_double_complex) , intent(out) :: den2p(1:nact, 1:nact, 1:nact, 1:nact, 0:(nproc-1))

  integer(c_long), external :: util_omp_iproc
  integer(c_long) :: istr,jstr,kstr,lstr,i1x_m,i1x,j1x_m,j1x,m_ij,m_kl,iact,jact,kact,lact, &
       & iord,jord,idist,jdist,kdist,ldist,distb,lla,ula,iproc,mvala,sfac,ssgn

  if (smul == 1 .and. nelact(1) == nelact(2)) then
!  if (.false.) then
     sfac = 2
  else
     sfac = 1
  end if

  !$omp parallel default(shared) private(mvala,idist,kdist,ldist,m_kl,i1x,iact,jact,kstr, &
  !$omp & ssgn,iord,lla,ula,j1x,kact,lact,lstr,jord,distb,iproc)
  iproc = util_omp_iproc()
  !$omp do
  do istr = 1, nstr_beta
     mvala = mtot-mval_beta(istr)
     if (mvala>mmax_alph .or. mvala<mmin_alph) cycle

     idist = dist_str_beta(1, istr)
     do m_ij = -mmax2, mmax2
        m_kl = -m_ij
        do i1x_m = 1, n1x_m_beta(m_ij,istr)
           i1x = map1x_m_beta(i1x_m,m_ij,istr)
           iact = h1x_beta (i1x,istr)
           jact = p1x_beta (i1x,istr)
           !20171201
           if (iact >= jact) cycle
           !20171201
           kstr = eq1x_beta(i1x,istr)
           kdist = dist_str_beta(1,kstr)
           ssgn = sgn1x_beta(i1x,istr) * sfac
!2RDMx           iord = nact * (iact - 1) + jact

           do jdist = 1, ndist_alph

!2RDM         if (det_allowed(jdist,idist) == 0) cycle
!2RDMx
              if (det_allowed(jdist,idist) == 0 .or. &
                  det_allowed(jdist,kdist) == 0) cycle
!2RDMx
              lla = llstr_dist_m_alph(jdist,mvala)
              ula =  nstr_dist_m_alph(jdist,mvala) + lla - 1
              do jstr = lla, ula
!                 work(jstr,iproc) = conjg(cic(mapf_detx(jstr,istr))) * ssgn
                 work(jstr,iproc) = conjg(cic(ntot_alph_beta(istr)+jstr)) * ssgn
              end do
           end do

           do j1x_m = 1, n1x_m_beta(m_kl,kstr)
              j1x = map1x_m_beta(j1x_m,m_kl,kstr)
              kact = h1x_beta (j1x,kstr)
              lact = p1x_beta (j1x,kstr)
              !20171201
              if (kact <= lact) cycle
              !20171201
              lstr = eq1x_beta(j1x,kstr)
              ldist = dist_str_beta(1,lstr)
!2RDMx              jord = nact * (kact - 1) + lact
              distb = max(max(idist,ldist),kdist)

!2RDMx              if (jord > iord) cycle
              do jdist = 1, ndist_alph
                 if (det_allowed(jdist,distb) == 0) cycle
                 lla = llstr_dist_m_alph(jdist,mvala)
                 ula =  nstr_dist_m_alph(jdist,mvala) + lla - 1
                 do jstr = lla, ula
                    den2p(jact, iact, kact, lact, iproc) = &
                    den2p(jact, iact, kact, lact, iproc) + &
!                    work(jstr, iproc) * sgn1x_beta(j1x, kstr) * cic(mapf_detx(jstr,lstr))
                    work(jstr, iproc) * sgn1x_beta(j1x, kstr) * cic(ntot_alph_beta(lstr)+jstr)
                 end do
              end do
           end do
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ormas_mkden2x_ras_bbp
!######################################################################
!######################################################################
subroutine ormas_mkden2x_ras_aap(cic, lwork, work, den2p, nproc)

  use, intrinsic :: iso_c_binding
  use mod_sph, only : mmax2
  use mod_bas, only : mtot, smul
  use mod_ormas, only : ndetx,ntot_alph_beta
  use mod_ormas, only : nact,det_allowed,nstr_beta,nstr_alph,nelact
  use mod_ormas, only : ndist_beta,nstr_dist_m_beta,llstr_dist_m_beta
  use mod_ormas, only : dist_str_alph,n1x_m_alph,map1x_m_alph,mval_alph
  use mod_ormas, only : p1x_alph,h1x_alph,eq1x_alph,sgn1x_alph,mmin_beta,mmax_beta

  implicit none
  integer(c_long), intent(in) :: lwork, nproc
  complex(c_double_complex) , intent(in) :: cic(1:ndetx)
  complex(c_double_complex) , intent(out) :: work(1:lwork, 0:(nproc-1))
  complex(c_double_complex) , intent(out) :: den2p(1:nact, 1:nact, 1:nact, 1:nact, 0:(nproc-1))

  integer(c_long), external :: util_omp_iproc
  integer(c_long) :: istr,jstr,kstr,lstr,i1x_m,i1x,j1x_m,j1x,m_ij,m_kl,iact,jact,kact,lact, &
       & iord,jord,idist,jdist,kdist,ldist,llb,ulb,iproc,mvalb,ssgn,dista

  if (smul == 1 .and. nelact(1) == nelact(2)) return
!  if (.false.) return

  !$omp parallel default(shared) private(mvalb,idist,kdist,ldist,m_kl,i1x,iact,jact,kstr, &
  !$omp & ssgn,iord,llb,ulb,j1x,kact,lact,lstr,jord,dista,iproc)
  iproc = util_omp_iproc()
  !$omp do
  do istr = 1, nstr_alph
     mvalb = mtot-mval_alph(istr)
     if (mvalb>mmax_beta .or. mvalb<mmin_beta) cycle

     idist = dist_str_alph(1, istr)
     do m_ij = -mmax2, mmax2
        m_kl = -m_ij
        do i1x_m = 1, n1x_m_alph(m_ij,istr)
           i1x = map1x_m_alph(i1x_m,m_ij,istr)
           iact = h1x_alph (i1x,istr)
           jact = p1x_alph (i1x,istr)
           !20171201
           if (iact >= jact) cycle
           !20171201
           kstr = eq1x_alph(i1x,istr)
           kdist = dist_str_alph(1,kstr)
           ssgn = sgn1x_alph(i1x,istr)
!2RDMx           iord = nact * (iact - 1) + jact

           do jdist = 1, ndist_beta

!2RDM         if (det_allowed(idist,jdist) == 0) cycle
!2RDMx
              if (det_allowed(idist,jdist) == 0 .or. &
                  det_allowed(kdist,jdist) == 0) cycle
!2RDMx

              llb = llstr_dist_m_beta(Jdist,mvalb)
              ulb =  nstr_dist_m_beta(jdist,mvalb) + llb - 1
              do jstr = llb, ulb
!                 work(jstr,iproc) = conjg(cic(mapf_detx(istr,jstr))) * ssgn
                 work(jstr,iproc) = conjg(cic(ntot_alph_beta(jstr)+istr)) * ssgn
              end do
           end do

           do j1x_m = 1, n1x_m_alph(m_kl,kstr)
              j1x = map1x_m_alph(j1x_m,m_kl,kstr)
              kact = h1x_alph (j1x,kstr)
              lact = p1x_alph (j1x,kstr)
              !20171201
              if (kact <= lact) cycle
              !20171201
              lstr = eq1x_alph(j1x,kstr)
              ldist = dist_str_alph(1,lstr)
!2RDMx              jord = nact * (kact - 1) + lact
              dista = max(max(idist,ldist),kdist)

!2RDMx              if (jord > iord) cycle
              do jdist = 1, ndist_beta
                 if (det_allowed(dista,jdist) == 0) cycle
                 llb = llstr_dist_m_beta(jdist,mvalb)
                 ulb =  nstr_dist_m_beta(jdist,mvalb) + llb - 1
                 do jstr = llb, ulb
                    den2p(jact, iact, kact, lact, iproc) = &
                    den2p(jact, iact, kact, lact, iproc) + &
!                    work(jstr, iproc) * sgn1x_alph(j1x, kstr) * cic(mapf_detx(lstr,jstr))
                    work(jstr, iproc) * sgn1x_alph(j1x, kstr) * cic(ntot_alph_beta(jstr)+lstr)
                 end do
              end do
           end do
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ormas_mkden2x_ras_aap
!######################################################################
!######################################################################
subroutine ormas_mkden2x_ras_abp(cic, lwork, work, den2p, nproc)

  use, intrinsic :: iso_c_binding
  use mod_sph, only : mmax2
  use mod_bas, only : mtot, smul
  use mod_const, only : czero,runit,ctwo
  use mod_ormas, only : ndetx,ntot_alph_beta
  use mod_ormas, only : nact,det_allowed,nstr_alph_beta,llstr_alph_beta,nelact
  use mod_ormas, only : ndist_alph,nstr_dist_m_alph,llstr_dist_m_alph,mmin_alph,mmax_alph
  use mod_ormas, only : ndist_beta,nstr_dist_m_beta,llstr_dist_m_beta,mmin_beta,mmax_beta
  use mod_ormas, only : dist_str_beta,n1x_m_beta,map1x_m_beta,mval_alph,mval_beta,dist_str_alph
  use mod_ormas, only : nstr_alph,n1x_m_alph,map1x_m_alph,p1x_alph,h1x_alph,eq1x_alph,sgn1x_alph
  use mod_ormas, only : nstr_beta,n1x_m_beta,map1x_m_beta,p1x_beta,h1x_beta,eq1x_beta,sgn1x_beta

  implicit none
  integer(c_long), intent(in) :: lwork, nproc
  complex(c_double_complex) , intent(in) :: cic(1:ndetx)
  complex(c_double_complex) , intent(out) :: work(1:lwork, 0:(nproc-1))
  complex(c_double_complex) , intent(out) :: den2p(1:nact, 1:nact, 1:nact, 1:nact, 0:(nproc-1))

  complex(c_double_complex) :: tmp
  complex(c_double_complex), allocatable :: d2cp(:,:,:,:)
  integer(c_long), external :: util_omp_iproc
  integer(c_long) :: istr,jstr,kstr,lstr,i1x_m,i1x,j1x_m,j1x,m_ij,m_kl,iact,jact,kact,lact, &
       & iord,jord,idist,jdist,kdist,ldist,distb,dista,lla,ula,llb,ulb,iproc,mvala,mvalb,sfac,ssgn

  !$omp parallel default(shared) private(mvala,idist,kdist,ldist,m_kl,i1x,iact,jact,kstr, &
  !$omp & ssgn,iord,lla,ula,j1x,kact,lact,lstr,jord,distb,tmp,iproc)
  iproc = util_omp_iproc()
  !$omp do
  do istr = 1, nstr_beta
     mvala = mtot-mval_beta(istr)
     if (mvala>mmax_alph .or. mvala<mmin_alph) cycle

     idist = dist_str_beta(1,istr)
     do m_ij = -mmax2, mmax2
        m_kl = -m_ij
        do i1x_m = 1, n1x_m_beta(m_ij,istr)
           i1x = map1x_m_beta(i1x_m,m_ij,istr)
           iact = h1x_beta(i1x,istr)
           jact = p1x_beta(i1x,istr)
           !20171201
           if (iact >= jact) cycle
           !20171201
           kstr = eq1x_beta(i1x,istr)
           kdist = dist_str_beta(1,kstr)
           ssgn = sgn1x_beta(i1x,istr)
!2RDMx           iord = nact * (iact - 1) + jact
!          distb = max(idist,kdist)

           do jdist = 1, ndist_alph
!             if (det_allowed(jdist,distb) == 0) cycle
              if (det_allowed(jdist,idist) == 0 .or. &
                  det_allowed(jdist,kdist) == 0) cycle
              lla = llstr_dist_m_alph(jdist,mvala)
              ula =  nstr_dist_m_alph(jdist,mvala) + lla - 1
              do jstr = lla, ula
!                 tmp = conjg(cic(mapf_detx(jstr,istr))) * ssgn
                 tmp = conjg(cic(ntot_alph_beta(istr)+jstr)) * ssgn
                 do j1x_m = 1, n1x_m_alph(m_kl,jstr)
                    j1x = map1x_m_alph(j1x_m,m_kl,jstr)
                    kact = h1x_alph(j1x,jstr)
                    lact = p1x_alph(j1x,jstr)
                    !20171201
                    if (kact <= lact) cycle
                    !20171201
                    lstr = eq1x_alph(j1x,jstr)
                    ldist = dist_str_alph(1,lstr)
!2RDMx                    jord = nact * (kact - 1) + lact
                    if (det_allowed(ldist,kdist) == 0) cycle
                    den2p(jact, iact, kact, lact, iproc) = &
                    den2p(jact, iact, kact, lact, iproc) + &
!                    tmp * sgn1x_alph(j1x,jstr) * cic(mapf_detx(lstr,kstr))
                    tmp * sgn1x_alph(j1x,jstr) * cic(ntot_alph_beta(kstr)+lstr)
                 end do
              end do
           end do
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel

!ver1  !$omp parallel default(shared) private(i1x,iact,jact,lstr,iord,tmp,m_kl,j1x,kact,lact,kstr,jord,iproc)
!ver1  iproc = util_omp_iproc()
!ver1  !$omp do
!ver1  do istr = 1, nstr_beta
!ver1     do jstr = llstr_alph_beta(istr), llstr_alph_beta(istr)+nstr_alph_beta(istr)-1
!ver1        do m_ij = -mmax2, mmax2
!ver1           do i1x_m = 1, n1x_m_beta(m_ij,istr)
!ver1              i1x = map1x_m_beta(i1x_m,m_ij,istr)
!ver1              iact = h1x_beta (i1x,istr)
!ver1              jact = p1x_beta (i1x,istr)
!ver1              lstr = eq1x_beta(i1x,istr)
!ver1              iord = nact * (iact - 1) + jact
!ver1              if (det_allowed(dist_str_alph(1, jstr),dist_str_beta(1, lstr)) == 0) cycle
!ver1
!ver1
!ver1              tmp = conjg(cic(mapf_detx(jstr,istr))) * sgn1x_beta(i1x,istr)
!ver1              m_kl = -m_ij
!ver1              do j1x_m = 1, n1x_m_alph(m_kl,jstr)
!ver1                 j1x = map1x_m_alph(j1x_m,m_kl,jstr)
!ver1                 kact = h1x_alph (j1x,jstr)
!ver1                 lact = p1x_alph (j1x,jstr)
!ver1                 kstr = eq1x_alph(j1x,jstr)
!ver1                 jord = nact * (kact - 1) + lact
!ver1!2RDMx                 if (jord > iord) cycle
!ver1
!ver1                 if (det_allowed(dist_str_alph(1, kstr),dist_str_beta(1, lstr)) == 0) cycle
!ver1        
!ver1
!ver1                 den2p(jact, iact, kact, lact, iproc) = &
!ver1                 den2p(jact, iact, kact, lact, iproc) + &
!ver1                 tmp * sgn1x_beta(j1x,jstr) * cic(mapf_detx(kstr,lstr))
!ver1              end do
!ver1           end do
!ver1        end do
!ver1     end do
!ver1  end do
!ver1  !$omp end do
!ver1  !$omp end parallel
  if (smul == 1 .and. nelact(1) == nelact(2)) then
!  if (.false.) then
     den2p = den2p + den2p
  else
     !$omp parallel default(shared) private(mvalb,idist,kdist,ldist,m_kl,i1x,iact,jact,kstr, &
     !$omp & ssgn,iord,llb,ulb,j1x,kact,lact,lstr,jord,dista,tmp,iproc)
     iproc = util_omp_iproc()
!20171201
     do iact = 1, nact
        do jact = 1, nact
           do kact = 1, nact
              do lact = 1, nact
                 den2p(jact, iact, kact, lact, iproc) = &
                 den2p(jact, iact, kact, lact, iproc) + &
           conjg(den2p(kact, lact, jact, iact, iproc))
              end do
           end do
        end do
     end do
!20171201     !$omp do
!20171201     do istr = 1, nstr_alph
!20171201        mvalb = mtot-mval_alph(istr)
!20171201        if (mvalb>mmax_beta .or. mvalb<mmin_beta) cycle
!20171201
!20171201        idist = dist_str_alph(1,istr)
!20171201        do m_ij = -mmax2, mmax2
!20171201           m_kl = -m_ij
!20171201           do i1x_m = 1, n1x_m_alph(m_ij,istr)
!20171201              i1x = map1x_m_alph(i1x_m,m_ij,istr)
!20171201              iact = h1x_alph(i1x,istr)
!20171201              jact = p1x_alph(i1x,istr)
!20171201              !20171201
!20171201              if (iact >= jact) cycle
!20171201              !20171201
!20171201              kstr = eq1x_alph(i1x,istr)
!20171201              kdist = dist_str_alph(1,kstr)
!20171201              ssgn = sgn1x_alph(i1x,istr)
!20171201!2RDMx              iord = nact * (iact - 1) + jact
!20171201!             dista = max(idist,kdist)
!20171201     
!20171201              do jdist = 1, ndist_beta
!20171201!                if (det_allowed(dista,jdist) == 0) cycle
!20171201                 if (det_allowed(idist,jdist) == 0 .or. &
!20171201                     det_allowed(kdist,jdist) == 0) cycle
!20171201                 llb = llstr_dist_m_beta(jdist,mvalb)
!20171201                 ulb =  nstr_dist_m_beta(jdist,mvalb) + llb - 1
!20171201                 do jstr = llb, ulb
!20171201!                    tmp = conjg(cic(mapf_detx(istr,jstr))) * ssgn
!20171201                    tmp = conjg(cic(ntot_alph_beta(jstr)+istr)) * ssgn
!20171201                    do j1x_m = 1, n1x_m_beta(m_kl,jstr)
!20171201                       j1x = map1x_m_beta(j1x_m,m_kl,jstr)
!20171201                       kact = h1x_beta(j1x,jstr)
!20171201                       lact = p1x_beta(j1x,jstr)
!20171201                       !20171201
!20171201                       if (kact <= lact) cycle
!20171201                       !20171201
!20171201                       lstr = eq1x_beta(j1x,jstr)
!20171201                       ldist = dist_str_beta(1,lstr)
!20171201!2RDMx                       jord = nact * (kact - 1) + lact
!20171201                       if (det_allowed(kdist,ldist) == 0) cycle
!20171201                       den2p(jact, iact, kact, lact, iproc) = &
!20171201                       den2p(jact, iact, kact, lact, iproc) + &
!20171201!                       tmp * sgn1x_beta(j1x,jstr) * cic(mapf_detx(kstr,lstr))
!20171201                       tmp * sgn1x_beta(j1x,jstr) * cic(ntot_alph_beta(lstr)+kstr)
!20171201                    end do
!20171201                 end do
!20171201              end do
!20171201           end do
!20171201        end do
!20171201     end do
!20171201     !$omp end do
     !$omp end parallel
!nyi     allocate(d2cp(1:nact,1:nact,1:nact,1:nact))
!nyi     do iproc = 0, nproc - 1
!nyi        d2cp(:,:,:,:) = den2p(:,:,:,:,iproc)
!nyi        do iact = 1, nact
!nyi        do jact = 1, nact
!nyi           do kact = 1, nact
!nyi           do lact = 1, nact
!nyi              jord = nact * (kact - 1) + lact
!nyi              if (jord > iord) cycle
!nyi              den2p(jact, iact, lact, kact, iproc) = &
!nyi              den2p(jact, iact, lact, kact, iproc) + &
!nyi              d2cp (lact, kact, jact, iact)
!nyi           end do
!nyi           end do
!nyi        end do
!nyi        end do
!nyi     end do
!nyi     deallocate(d2cp)
  end if

end subroutine ormas_mkden2x_ras_abp
!######################################################################
!######################################################################
subroutine ormas_mkden2x_ras_sum(den2aa,den2bb,den2ab,den2,nproc)

  use, intrinsic :: iso_c_binding
  use mod_const, only : runit
  use mod_ormas, only : nact,ncore
  use mod_bas, only : mval

  implicit none
  integer(c_long), intent(in) :: nproc
  complex(c_double_complex) , intent(in) :: den2aa(1:nact, 1:nact, 1:nact, 1:nact, 0:(nproc-1))
  complex(c_double_complex) , intent(in) :: den2bb(1:nact, 1:nact, 1:nact, 1:nact, 0:(nproc-1))
  complex(c_double_complex) , intent(in) :: den2ab(1:nact, 1:nact, 1:nact, 1:nact, 0:(nproc-1))
  complex(c_double_complex) , intent(out) :: den2(1:nact, 1:nact, 1:nact, 1:nact)

  integer(c_long) :: iact, jact, kact, lact, iord, jord, iproc

  do iproc = 0, nproc - 1
     den2(:,:,:,:) = den2(:,:,:,:) + den2aa(:,:,:,:,iproc) &
                                   + den2bb(:,:,:,:,iproc) &
                                   + den2ab(:,:,:,:,iproc)
  end do

  !20171201
  do iact = 1, nact
  do jact = 1, nact
  do kact = 1, nact
  do lact = 1, nact
     if (mval(ncore+iact).ne.mval(ncore+jact).or. &
         mval(ncore+iact).ne.mval(ncore+jact)) then
        den2(iact,jact,kact,lact) = 0d0
     end if
  end do
  end do
  end do
  end do
  !20171201

!2RDMx  ! particle inter-change symmetry
!2RDMx  do iact = 1, nact
!2RDMx  do jact = 1, nact
!2RDMx     iord = nact * (iact - 1) + jact
!2RDMx     do kact = 1, nact
!2RDMx     do lact = 1, nact
!2RDMx        jord = nact * (kact - 1) + lact
!2RDMx        if (jord <= iord) cycle
!2RDMx        den2(jact, iact, lact, kact) = &
!2RDMx        den2(lact, kact, jact, iact)
!2RDMx     end do
!2RDMx     end do
!2RDMx  end do
!2RDMx  end do

end subroutine ormas_mkden2x_ras_sum
!######################################################################
