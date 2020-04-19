!######################################################################
subroutine ormas_mkden2_ras(cic, den1, den2)

  use, intrinsic :: iso_c_binding
  use mod_const, only : czero
  use mod_ormas, only : nact,nstr_alph,nstr_beta,ndetx,nelact

  implicit none
  complex(c_double_complex) , intent(in) :: cic(1:ndetx)
  complex(c_double_complex) , intent(in) :: den1(1:nact, 1:nact)
  complex(c_double_complex) , intent(out) :: den2(1:nact, 1:nact, 1:nact, 1:nact)

  integer(c_long) :: lwork
  integer(c_long) :: nproc,iproc,iact,jact,kact,lact
  complex(c_double_complex) , allocatable :: work(:,:)
  complex(c_double_complex) , allocatable :: den2aa(:,:,:,:,:)
  complex(c_double_complex) , allocatable :: den2bb(:,:,:,:,:)
  complex(c_double_complex) , allocatable :: den2ab(:,:,:,:,:)
  integer(c_long), external :: util_omp_nproc

  if (nact==0) return

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

!bug!  if (nelact(2) >= 2) call ormas_mkden2_ras_bbp(cic,lwork,work,den2bb,nproc)
!bug!  if (nelact(1) >= 2) call ormas_mkden2_ras_aap(cic,lwork,work,den2aa,nproc)
  if (nelact(2) >= 1) call ormas_mkden2_ras_bbp(cic,lwork,work,den2bb,nproc)
  if (nelact(1) >= 1) call ormas_mkden2_ras_aap(cic,lwork,work,den2aa,nproc)
  if (nelact(1) >= 1 .and. nelact(2) >= 1) call ormas_mkden2_ras_abp(cic,lwork,work,den2ab,nproc)
  if (nelact(3) >= 2) call ormas_mkden2_ras_sum(den1,den2aa,den2bb,den2ab,den2,nproc)

  deallocate(den2ab)
  deallocate(den2bb)
  deallocate(den2aa)
  deallocate(work)

!DEBUG
!  call util_print_vec(nact**4, den2, "test.den2")
!DEBUG

end subroutine ormas_mkden2_ras
!######################################################################
!######################################################################
subroutine ormas_mkden2_ras_bbp(cic, lwork, work, den2p, nproc)

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
       & iord,jord,idist,jdist,ldist,lla,ula,iproc,mvala,sfac,ssgn

  if (smul == 1 .and. nelact(1) == nelact(2)) then
     sfac = 2
  else
     sfac = 1
  end if

  !$omp parallel default(shared) private(mvala,idist,ldist,m_kl,i1x,iact,jact,kstr, &
  !$omp & ssgn,iord,lla,ula,j1x,kact,lact,lstr,jord,iproc)
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
           iact = h1x_beta (i1x,istr)
           jact = p1x_beta (i1x,istr)
           kstr = eq1x_beta(i1x,istr)
           ssgn = sgn1x_beta(i1x,istr) * sfac
           iord = nact * (iact - 1) + jact

           do jdist = 1, ndist_alph
              if (det_allowed(jdist,idist) == 0) cycle
              lla = llstr_dist_m_alph(jdist,mvala)
              ula =  nstr_dist_m_alph(jdist,mvala) + lla - 1
              do jstr = lla, ula
!                work(jstr,iproc) = conjg(cic(mapf_detx(jstr,istr))) * ssgn
                 work(jstr,iproc) = conjg(cic(ntot_alph_beta(istr)+jstr)) * ssgn
              end do
           end do

           do j1x_m = 1, n1x_m_beta(m_kl,kstr)
              j1x = map1x_m_beta(j1x_m,m_kl,kstr)
              kact = h1x_beta (j1x,kstr)
              lact = p1x_beta (j1x,kstr)
              lstr = eq1x_beta(j1x,kstr)
              jord = nact * (kact - 1) + lact
              ldist = dist_str_beta(1,lstr)

              if (jord > iord) cycle
              do jdist = 1, ndist_alph
                 if (det_allowed(jdist,idist) == 0 .or. &
                     det_allowed(jdist,ldist) == 0) cycle
                 lla = llstr_dist_m_alph(jdist,mvala)
                 ula =  nstr_dist_m_alph(jdist,mvala) + lla - 1
                 do jstr = lla, ula
                    den2p(jact,iact,lact,kact,iproc) = &
                    den2p(jact,iact,lact,kact,iproc) + &
!                   work(jstr,iproc) * sgn1x_beta(j1x,kstr) * cic(mapf_detx(jstr,lstr))
                    work(jstr,iproc) * sgn1x_beta(j1x,kstr) * cic(ntot_alph_beta(lstr)+jstr)
                 end do
              end do
           end do
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ormas_mkden2_ras_bbp
!######################################################################
!######################################################################
subroutine ormas_mkden2_ras_aap(cic, lwork, work, den2p, nproc)

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
       & iord,jord,idist,jdist,ldist,llb,ulb,iproc,mvalb,ssgn

  if (smul == 1 .and. nelact(1) == nelact(2)) return
!  if (.false.) return

  !$omp parallel default(shared) private(mvalb,idist,ldist,m_kl,i1x,iact,jact,kstr, &
  !$omp & ssgn,iord,llb,ulb,j1x,kact,lact,lstr,jord,iproc)
  iproc = util_omp_iproc()
  !$omp do
  do istr = 1, nstr_alph
     mvalb = mtot-mval_alph(istr)
     if (mvalb>mmax_beta .or. mvalb<mmin_beta) cycle

     idist = dist_str_alph(1,istr)
     do m_ij = -mmax2, mmax2
        m_kl = -m_ij
        do i1x_m = 1, n1x_m_alph(m_ij,istr)
           i1x = map1x_m_alph(i1x_m,m_ij,istr)
           iact = h1x_alph (i1x,istr)
           jact = p1x_alph (i1x,istr)
           kstr = eq1x_alph(i1x,istr)
           ssgn = sgn1x_alph(i1x,istr)
           iord = nact * (iact - 1) + jact

           do jdist = 1, ndist_beta
              if (det_allowed(idist,jdist) == 0) cycle
              llb = llstr_dist_m_beta(jdist,mvalb)
              ulb =  nstr_dist_m_beta(jdist,mvalb) + llb - 1
              do jstr = llb, ulb
!                work(jstr,iproc) = conjg(cic(mapf_detx(istr,jstr))) * ssgn
                 work(jstr,iproc) = conjg(cic(ntot_alph_beta(jstr)+istr)) * ssgn
              end do
           end do

           do j1x_m = 1, n1x_m_alph(m_kl,kstr)
              j1x = map1x_m_alph(j1x_m,m_kl,kstr)
              kact = h1x_alph (j1x,kstr)
              lact = p1x_alph (j1x,kstr)
              lstr = eq1x_alph(j1x,kstr)
              jord = nact * (kact - 1) + lact
              ldist = dist_str_alph(1,lstr)

              if (jord > iord) cycle
              do jdist = 1, ndist_beta
                 if (det_allowed(idist,jdist) == 0 .or. &
                     det_allowed(ldist,jdist) == 0) cycle
                 llb = llstr_dist_m_beta(jdist,mvalb)
                 ulb =  nstr_dist_m_beta(jdist,mvalb) + llb - 1
                 do jstr = llb, ulb
                    den2p(jact,iact,lact,kact,iproc) = &
                    den2p(jact,iact,lact,kact,iproc) + &
!                   work(jstr,iproc) * sgn1x_alph(j1x,kstr) * cic(mapf_detx(lstr,jstr))
                    work(jstr,iproc) * sgn1x_alph(j1x,kstr) * cic(ntot_alph_beta(jstr)+lstr)
                 end do
              end do
           end do
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ormas_mkden2_ras_aap
!######################################################################
!######################################################################
subroutine ormas_mkden2_ras_abp(cic, lwork, work, den2p, nproc)

  use, intrinsic :: iso_c_binding
  use mod_sph, only : mmax2
  use mod_bas, only : smul,mtot
  use mod_const, only : czero,runit,ctwo
  use mod_ormas, only : ndetx,ntot_alph_beta
  use mod_ormas, only : dist_str_alph,dist_str_beta,mval_alph,mval_beta
  use mod_ormas, only : nact,det_allowed,nstr_alph_beta,llstr_alph_beta,nelact
  use mod_ormas, only : nstr_beta_alph,llstr_beta_alph,mmin_alph,mmax_alph,mmin_beta,mmax_beta
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
  integer(c_long) :: m_ij,m_kl,kdist,ldist,mvala,mvalb
  integer(c_long) :: istr,jstr,kstr,lstr,i1x_m,i1x,j1x_m,j1x,iact,jact,kact,lact,iord,jord,iproc

  !$omp parallel default(shared) private(mvala,i1x,iact,jact,lstr,iord, &
  !$omp & tmp,m_kl,j1x,kact,lact,kstr,jord,kdist,ldist,iproc)
  iproc = util_omp_iproc()
  !$omp do
  do istr = 1, nstr_beta
     mvala = mtot-mval_beta(istr)
     if (mvala>mmax_alph .or. mvala<mmin_alph) cycle
     do jstr = llstr_alph_beta(istr), llstr_alph_beta(istr)+nstr_alph_beta(istr)-1
        do m_ij = -mmax2, mmax2
           do i1x_m = 1, n1x_m_beta(m_ij,istr)
              i1x = map1x_m_beta(i1x_m,m_ij,istr)
              iact = h1x_beta (i1x,istr)
              jact = p1x_beta (i1x,istr)
              lstr = eq1x_beta(i1x,istr)
              ldist = dist_str_beta(1,lstr)
              iord = nact * (iact - 1) + jact
!             tmp = conjg(cic(mapf_detx(jstr,istr))) * sgn1x_beta(i1x,istr)
              tmp = conjg(cic(ntot_alph_beta(istr)+jstr)) * sgn1x_beta(i1x,istr)
              m_kl = -m_ij
              do j1x_m = 1, n1x_m_alph(m_kl,jstr)
                 j1x = map1x_m_alph(j1x_m,m_kl,jstr)
                 kact = h1x_alph (j1x,jstr)
                 lact = p1x_alph (j1x,jstr)
                 kstr = eq1x_alph(j1x,jstr)
                 kdist = dist_str_alph(1,kstr)
                 jord = nact * (kact - 1) + lact
!bug             if (jord > iord .or. det_allowed(kdist,ldist) == 0) cycle
                 if (det_allowed(kdist,ldist) == 0) cycle
                 den2p(jact,iact,lact,kact,iproc) = &
                 den2p(jact,iact,lact,kact,iproc) + &
!                tmp * sgn1x_alph(j1x,jstr) * cic(mapf_detx(kstr,lstr))
                 tmp * sgn1x_alph(j1x,jstr) * cic(ntot_alph_beta(lstr)+kstr)
              end do
           end do
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel

  if (smul == 1 .and. nelact(1) == nelact(2)) then
     den2p = den2p + den2p
  else
     allocate(d2cp(1:nact,1:nact,1:nact,1:nact))
     do iproc = 0, nproc - 1
        d2cp(:,:,:,:) = den2p(:,:,:,:,iproc)
        do iact = 1, nact
        do jact = 1, nact
           do kact = 1, nact
           do lact = 1, nact
              den2p(jact, iact, lact, kact, iproc) = &
              den2p(jact, iact, lact, kact, iproc) + &
              d2cp (lact, kact, jact, iact)
           end do
           end do
        end do
        end do
     end do
     deallocate(d2cp)
  end if

end subroutine ormas_mkden2_ras_abp
!######################################################################
!######################################################################
subroutine ormas_mkden2_ras_sum(den1,den2aa,den2bb,den2ab,den2,nproc)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : smul
  use mod_const, only : runit
  use mod_ormas, only : nact,nelact

  implicit none
  integer(c_long), intent(in) :: nproc
  complex(c_double_complex) , intent(in) :: den1(1:nact, 1:nact)
  complex(c_double_complex) , intent(inout) :: den2aa(1:nact, 1:nact, 1:nact, 1:nact, 0:(nproc-1))
  complex(c_double_complex) , intent(inout) :: den2bb(1:nact, 1:nact, 1:nact, 1:nact, 0:(nproc-1))
  complex(c_double_complex) , intent(inout) :: den2ab(1:nact, 1:nact, 1:nact, 1:nact, 0:(nproc-1))
  complex(c_double_complex) , intent(out) :: den2(1:nact, 1:nact, 1:nact, 1:nact)

  integer(c_long) :: iact, jact, kact, lact, iord, jord, iproc
  complex(c_double_complex) :: daa,dbb,dab

  do iproc = 0, nproc - 1
     den2(:,:,:,:) = den2(:,:,:,:) + den2aa(:,:,:,:,iproc) &
                                   + den2bb(:,:,:,:,iproc) &
                                   + den2ab(:,:,:,:,iproc)
  end do

  ! 1-rdm contributions
  do iact = 1, nact
     do jact = 1, nact
        do kact = 1, nact
           den2(jact, iact, kact, jact) = &
         & den2(jact, iact, kact, jact) - den1(kact, iact)
        end do
     end do
  end do

  ! particle inter-change symmetry
  do iact = 1, nact
  do jact = 1, nact
     iord = nact * (iact - 1) + jact
     do kact = 1, nact
     do lact = 1, nact
        jord = nact * (kact - 1) + lact
        if (jord <= iord) cycle
        den2(jact, iact, lact, kact) = &
        den2(lact, kact, jact, iact)
     end do
     end do
  end do
  end do

  !DEBUG
  !do iproc = 1, nproc - 1
  !   den2aa(:,:,:,:,0) = den2aa(:,:,:,:,0) + den2aa(:,:,:,:,iproc)
  !   den2bb(:,:,:,:,0) = den2bb(:,:,:,:,0) + den2bb(:,:,:,:,iproc)
  !   den2ab(:,:,:,:,0) = den2ab(:,:,:,:,0) + den2ab(:,:,:,:,iproc)
  !end do
  !if (smul == 1 .and. nelact(1) == nelact(2)) then
  !   den2bb(:,:,:,:,0) = den2bb(:,:,:,:,0)*0.5d0
  !   den2aa(:,:,:,:,0) = den2bb(:,:,:,:,0)
  !end if
  !
  !do iact = 1, nact
  !do jact = 1, nact
  !   iord = nact * (iact - 1) + jact
  !   do kact = 1, nact
  !   do lact = 1, nact
  !      jord = nact * (kact - 1) + lact
  !      if (jord <= iord) cycle
  !      den2aa(jact,iact,lact,kact,0) = den2aa(lact,kact,jact,iact,0)
  !      den2bb(jact,iact,lact,kact,0) = den2bb(lact,kact,jact,iact,0)
  !      den2ab(jact,iact,lact,kact,0) = den2ab(lact,kact,jact,iact,0)
  !   end do
  !   end do
  !end do
  !end do
  !
  !write(6,"('den2ab:')")
  !do iact = 1, nact
  !do jact = 1, nact
  !   do kact = 1, nact
  !   do lact = 1, nact
  !      write(6,"(4i5,f20.10)") jact,iact,lact,kact,dble(den2ab(jact,iact,lact,kact,0))
  !   end do
  !   end do
  !end do
  !end do
  !
  !daa = 0d0
  !dbb = 0d0
  !dab = 0d0
  !do iact = 1, nact
  !   do jact = 1, nact
  !      daa = daa + den2aa(iact,iact,jact,jact,0)
  !      dbb = dbb + den2bb(iact,iact,jact,jact,0)
  !      dab = dab + den2ab(iact,iact,jact,jact,0)
  !   end do
  !end do
  !write(6, "('ormas_mkden2_ras_sum: ', 3f20.10)") dble(daa)-nelact(1),dble(dbb)-nelact(2),dble(dab)
  !DEBUG


end subroutine ormas_mkden2_ras_sum
!######################################################################
