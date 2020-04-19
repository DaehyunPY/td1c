!######################################################################
subroutine ormas_hcic_ras(int1e, int2e, cic, hcic)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : smul
  use mod_ormas, only : nact,nelact,ndetx,ormas_donly

  implicit none
  complex(c_double_complex), intent(in) :: int1e(1:nact, 1:nact)
  complex(c_double_complex), intent(in) :: int2e(1:nact, 1:nact, 1:nact, 1:nact)
  complex(c_double_complex), intent(in) :: cic(1:ndetx)
  complex(c_double_complex), intent(inout) :: hcic(1:ndetx)
  integer(c_long) :: istr,jstr
  complex(c_double_complex), allocatable :: eff1e(:,:), eff2e(:,:,:,:)

  if (nact == 0) return
  allocate(eff1e(1:nact, 1:nact))
  allocate(eff2e(1:nact, 1:nact, 1:nact, 1:nact))

  ! 1e integrals with separable 2e contributions
  call ormas_hcic_ras_effint1e(int1e, int2e, eff1e)
  call ormas_hcic_ras_effint2e(int2e, eff2e)

  if (nelact(2) >= 1) call ormas_hcic_ras_bbp(eff1e, eff2e, cic, hcic)
  if (nelact(1) >= 1) call ormas_hcic_ras_aap(eff1e, eff2e, cic, hcic)
  if (nelact(1) >= 1 .and. nelact(2) >= 1) call ormas_hcic_ras_abp(int2e, cic, hcic)

  deallocate(eff2e)
  deallocate(eff1e)

  if (ormas_donly) call ormas_rm_singles(hcic)

!DEBUG
!  call ormas_cic_print(hcic, "test.hcic")
!DEBUG

end subroutine ormas_hcic_ras
!######################################################################
subroutine ormas_hcic_ras_bbp(int1e, int2e, cic, hcic)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : thrcic
  use mod_const, only : czero
  use mod_sph, only : mmax2
  use mod_bas, only : mtot
  use mod_ormas, only : nact, det_allowed, nelact, ndetx, ntot_alph_beta
  use mod_ormas, only : ndist_alph, nstr_alph, lstr_alph_dist
  use mod_ormas, only : ndist_beta, nstr_beta, lstr_beta_dist, &
       & dist_str_beta, n1x_beta, p1x_beta, h1x_beta, eq1x_beta, sgn1x_beta
  use mod_ormas, only : n1x_m_alph, map1x_m_alph
  use mod_ormas, only : n1x_m_beta, map1x_m_beta
  use mod_ormas, only : nstr_dist_m_alph, llstr_dist_m_alph, mval_alph
  use mod_ormas, only : nstr_dist_m_beta, llstr_dist_m_beta, mval_beta,mmin_alph,mmax_alph

  implicit none
  complex(c_double_complex), intent(in) :: int1e(1:nact, 1:nact)
  complex(c_double_complex), intent(in) :: int2e(1:nact, 1:nact, 1:nact, 1:nact)
  complex(c_double_complex), intent(in) :: cic(1:ndetx)
  complex(c_double_complex), intent(inout) :: hcic(1:ndetx)

  complex(c_double_complex), allocatable :: sint(:,:)
  integer(c_long) :: istr, jstr, kstr, lstr, i1x, j1x, ifun, jfun, kfun, lfun, &
       & iord, jord, idist, jdist, kdist, ii, lla, ula, llb, ulb, iproc, nproc, &
       & i1x_m, j1x_m, m_ij, m_kl, tsgn, mvala
  integer(c_long), external :: util_omp_nproc
  integer(c_long), external :: util_omp_iproc

  nproc = util_omp_nproc()
  allocate(sint(1:nstr_beta, 0:(nproc-1)))

  !$omp parallel default(shared) &
  !$omp private(iproc,idist,kstr,lstr,ifun,jfun,kfun,lfun, &
  !$omp & iord,jord,lla,ula,llb,ulb,tsgn,i1x,j1x,m_ij,m_kl,mvala)
  iproc = util_omp_iproc()
  !$omp do
  do istr = 1, nstr_beta
     mvala = mtot-mval_beta(istr)
     if (mvala>mmax_alph .or. mvala<mmin_alph) cycle

     idist = dist_str_beta(1, istr)
     call util_zcopy(nstr_beta, czero, 0, sint(1,iproc), 1)
     do m_ij = -mmax2, mmax2
        m_kl = -m_ij
        do i1x_m = 1, n1x_m_beta(m_ij, istr)
           i1x = map1x_m_beta(i1x_m, m_ij, istr)
           ifun = h1x_beta  (i1x, istr)
           jfun = p1x_beta  (i1x, istr)
           kstr = eq1x_beta (i1x, istr)
           tsgn = sgn1x_beta(i1x, istr)
           iord = nact * (ifun - 1) + jfun
           if (m_ij == 0) sint(kstr,iproc) = sint(kstr,iproc) + tsgn * int1e(ifun,jfun)
           do j1x_m = 1, n1x_m_beta(m_kl,kstr)
              j1x = map1x_m_beta(j1x_m,m_kl,kstr)
              kfun = h1x_beta (j1x,kstr)
              lfun = p1x_beta (j1x,kstr)
              lstr = eq1x_beta(j1x,kstr)
              jord = nact * (kfun - 1) + lfun
              if (jord <= iord) then
                 sint(lstr,iproc) = sint(lstr,iproc) &
                      + tsgn * sgn1x_beta(j1x,kstr) * int2e(ifun,jfun,kfun,lfun)
              end if
           end do
        end do
     end do

     do kdist = 1, ndist_beta
        llb = llstr_dist_m_beta(kdist,mval_beta(istr))
        ulb =  nstr_dist_m_beta(kdist,mval_beta(istr)) + llb - 1
        do jdist = 1, ndist_alph
           if (det_allowed(jdist,kdist) == 0 .or. &
               det_allowed(jdist,idist) == 0) cycle
           lla = llstr_dist_m_alph(jdist,mvala)
           ula =  nstr_dist_m_alph(jdist,mvala) + lla - 1
           do jstr = llb, ulb
              if (abs(sint(jstr,iproc)) > thrcic) then
                 do ii = lla, ula
!1                    hcic(ii,istr) = hcic(ii,istr) + sint(jstr,iproc) * cic(ii,jstr)
!2                    hcic(mapf_detx(ii,istr)) = &
!2                    hcic(mapf_detx(ii,istr)) + sint(jstr,iproc) * &
!2                     cic(mapf_detx(ii,jstr))
                    hcic(ntot_alph_beta(istr)+ii) = &
                    hcic(ntot_alph_beta(istr)+ii) + sint(jstr,iproc) * &
                     cic(ntot_alph_beta(jstr)+ii)
                 end do
              end if
           end do
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel

  deallocate(sint)

end subroutine ormas_hcic_ras_bbp
!######################################################################
!######################################################################
subroutine ormas_hcic_ras_aap(int1e, int2e, cic, hcic)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : thrcic
  use mod_const, only : czero
  use mod_bas, only : smul, mtot
  use mod_sph, only : mmax2
  use mod_ormas, only : nact, det_allowed, nelact, ndetx, mapr_detx, ntot_alph_beta
  use mod_ormas, only : ndist_alph, nstr_alph
  use mod_ormas, only : ndist_beta, nstr_beta
  use mod_ormas, only : dist_str_alph, p1x_alph, h1x_alph, eq1x_alph, sgn1x_alph
  use mod_ormas, only : n1x_m_alph, map1x_m_alph,mmin_beta,mmax_beta
  use mod_ormas, only : nstr_dist_m_alph, llstr_dist_m_alph, mval_alph
  use mod_ormas, only : nstr_dist_m_beta, llstr_dist_m_beta, mval_beta

  implicit none
  complex(c_double_complex), intent(in) :: int1e(1:nact, 1:nact)
  complex(c_double_complex), intent(in) :: int2e(1:nact, 1:nact, 1:nact, 1:nact)
  complex(c_double_complex), intent(in) :: cic(1:ndetx)
  complex(c_double_complex), intent(inout) :: hcic(1:ndetx)

  complex(c_double_complex), allocatable :: sint(:,:), dcic(:)
  integer(c_long) :: istr, jstr, kstr, lstr, i1x, j1x, ifun, jfun, kfun, lfun, &
       & iord, jord, idist, jdist, kdist, ii, lla, ula, llb, ulb, iproc, nproc, &
       & i1x_m, j1x_m, m_ij, m_kl, tsgn, idet, jdet, mvalb
  integer(c_long), external :: util_omp_nproc
  integer(c_long), external :: util_omp_iproc

  !##### Applying the following up-down spin symmetry #####
  !##### leads to LESS stable propagation ! Why ??    #####
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
     !$omp end do
!!     !$omp do   
!!     !$omp parallel default(shared) private(istr,jstr)
!!     !$omp do
!!     do istr = 1, nstr_beta
!!        do jstr = 1, istr
!!           hcic(jstr,istr) = hcic(jstr,istr) + hcic(istr,jstr)
!!        end do
!!     end do
!!     !$omp end do
!!     !$omp do
!!     do istr = 1, nstr_beta
!!        do jstr = istr + 1, nstr_alph
!!           hcic(jstr,istr) = hcic(istr,jstr)
!!        end do
!!     end do
!!     !$omp end do
     !$omp end parallel
     deallocate(dcic)
     return
  end if

  nproc = util_omp_nproc()
  allocate(sint(1:nstr_alph, 0:(nproc-1)))
  !$omp parallel default(shared) &
  !$omp private(iproc,idist,kstr,lstr,ifun,jfun,kfun,lfun, &
  !$omp & iord,jord,lla,ula,llb,ulb,tsgn,i1x,j1x,m_ij,m_kl,mvalb)
  iproc = util_omp_iproc()
  !$omp do
  do istr = 1, nstr_alph
     mvalb = mtot-mval_alph(istr)
     if (mvalb>mmax_beta .or. mvalb<mmin_beta) cycle

     idist = dist_str_alph(1, istr)
     call util_zcopy(nstr_alph, czero, 0, sint(1,iproc), 1)
     do m_ij = -mmax2, mmax2
        m_kl = -m_ij
        do i1x_m = 1, n1x_m_alph(m_ij, istr)
           i1x = map1x_m_alph(i1x_m, m_ij, istr)
           ifun = h1x_alph  (i1x, istr)
           jfun = p1x_alph  (i1x, istr)
           kstr = eq1x_alph (i1x, istr)
           tsgn = sgn1x_alph(i1x, istr)
           iord = nact * (ifun - 1) + jfun
           if (m_ij == 0) sint(kstr, iproc) = sint(kstr, iproc) &
                + tsgn * int1e(ifun, jfun)
           do j1x_m = 1, n1x_m_alph(m_kl, kstr)
              j1x = map1x_m_alph(j1x_m, m_kl, kstr)
              kfun = h1x_alph (j1x, kstr)
              lfun = p1x_alph (j1x, kstr)
              lstr = eq1x_alph(j1x, kstr)
              jord = nact * (kfun - 1) + lfun
              if (jord <= iord) then
                 sint(lstr, iproc) = sint(lstr, iproc) &
                      + tsgn * sgn1x_alph(j1x, kstr) * int2e(ifun, jfun, kfun, lfun)
              end if
           end do
        end do
     end do

     do kdist = 1, ndist_alph
        lla = llstr_dist_m_alph(kdist,mval_alph(istr))
        ula =  nstr_dist_m_alph(kdist,mval_alph(istr)) + lla - 1
        do jdist = 1, ndist_beta
           if (det_allowed(kdist,jdist) == 0 .or. &
               det_allowed(idist,jdist) == 0) cycle
           llb = llstr_dist_m_beta(jdist, mvalb)
           ulb =  nstr_dist_m_beta(jdist, mvalb) + llb - 1
           do jstr = lla, ula
              if (abs(sint(jstr, iproc)) > thrcic) then
                 do ii = llb, ulb
!1                    hcic(istr,ii) = hcic(istr,ii) + sint(jstr,iproc) * cic(jstr,ii)
!2                    hcic(mapf_detx(istr,ii)) = &
!2                    hcic(mapf_detx(istr,ii)) + sint(jstr,iproc) * &
!2                     cic(mapf_detx(jstr,ii))
                    hcic(ntot_alph_beta(ii)+istr) = &
                    hcic(ntot_alph_beta(ii)+istr) + sint(jstr,iproc) * &
                     cic(ntot_alph_beta(ii)+jstr)
                 end do
              end if
           end do
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel
  deallocate(sint)

end subroutine ormas_hcic_ras_aap
!######################################################################
!######################################################################
subroutine ormas_hcic_ras_abp(int2e, cic, hcic)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : mtot
  use mod_sph, only : mmax2
  use mod_const, only : czero
  use mod_ormas, only : nact, det_allowed, ndist_alph, ndist_beta, ndetx, mapr_detx
  use mod_ormas, only : nstr_alph, n1x_alph, p1x_alph, h1x_alph, eq1x_alph, sgn1x_alph
  use mod_ormas, only : nstr_beta, n1x_beta, p1x_beta, h1x_beta, eq1x_beta, sgn1x_beta
  use mod_ormas, only : n1x_m_alph, map1x_m_alph, llstr_alph_beta, nstr_alph_beta
  use mod_ormas, only : n1x_m_beta, map1x_m_beta, dist_str_alph, dist_str_beta
  use mod_ormas, only : nstr_dist_m_alph, llstr_dist_m_alph, mval_alph
  use mod_ormas, only : nstr_dist_m_beta, llstr_dist_m_beta, mval_beta
  use mod_ormas, only : lstr_alph_dist, ntot_alph_beta, ntot_beta_alph

  implicit none
  complex(c_double_complex), intent(in) :: int2e(1:nact, 1:nact, 1:nact, 1:nact)
  complex(c_double_complex), intent(in) :: cic(1:ndetx)
  complex(c_double_complex), intent(inout) :: hcic(1:ndetx)
  integer(c_long) :: istr,jstr,kstr,lstr,ifun,jfun,kfun,lfun,i1x_m,i1x,j1x_m,j1x,m_ij,m_kl
  integer(c_long) :: idist,jdist,kdist,ldist,lla,ula,idet
  complex(c_double_complex) :: tmp

  !$omp parallel default(shared) private(jstr,istr,m_kl,i1x,ifun,jfun,j1x,kfun,lfun,kstr,lstr,tmp)
  !$omp do
  do idet = 1, ndetx
     jstr = mapr_detx(1,idet)
     istr = mapr_detx(2,idet)
     if (mval_alph(jstr)+mval_beta(istr).ne.mtot) cycle
     !##### DEBUG #####
     if (det_allowed(dist_str_alph(1,jstr),dist_str_beta(1,istr)) == 0) cycle
     !##### DEBUG #####
     do m_ij = -mmax2, mmax2
        do i1x_m = 1, n1x_m_beta(m_ij,istr)
           i1x = map1x_m_beta(i1x_m,m_ij,istr)
           ifun = h1x_beta  (i1x,istr)
           jfun = p1x_beta  (i1x,istr)
           lstr = eq1x_beta (i1x,istr)
           tmp = czero
           m_kl = -m_ij
           do j1x_m = 1, n1x_m_alph(m_kl, jstr)
              j1x = map1x_m_alph(j1x_m, m_kl, jstr)
              kfun = h1x_alph (j1x, jstr)
              lfun = p1x_alph (j1x, jstr)
              kstr = eq1x_alph(j1x, jstr)
              !### intensive if statement is not appropriate #####
              !### n1x_dist_m_aplh/beta should be developed ######
              if (det_allowed(dist_str_alph(1,kstr),dist_str_beta(1,lstr)) == 0) cycle
              !################################################
!2              tmp = tmp + sgn1x_alph(j1x,jstr)*int2e(ifun,jfun,kfun,lfun)*cic(mapf_detx(kstr,lstr))
              tmp = tmp + sgn1x_alph(j1x,jstr)*int2e(ifun,jfun,kfun,lfun)*cic(ntot_alph_beta(lstr)+kstr)
           end do
           hcic(idet) = hcic(idet) + tmp * sgn1x_beta(i1x,istr)
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel


!v4  !$omp parallel default(shared) private(idist,m_kl,i1x,ifun,jfun,&
!v4  !$omp & lstr,ldist,lla,ula,j1x,kfun,lfun,kstr,kdist,tmp)
!v4  !$omp do
!v4  do istr = 1, nstr_beta
!v4     idist = dist_str_beta(1,istr)
!v4     do m_ij = -mmax2, mmax2
!v4        m_kl = -m_ij
!v4        do i1x_m = 1, n1x_m_beta(m_ij,istr)
!v4           i1x = map1x_m_beta(i1x_m,m_ij,istr)
!v4           ifun = h1x_beta  (i1x,istr)
!v4           jfun = p1x_beta  (i1x,istr)
!v4           lstr = eq1x_beta (i1x,istr)
!v4           ldist = dist_str_beta(1,lstr)
!v4           do jdist = 1, ndist_alph
!v4              if (det_allowed(jdist,idist) == 0) cycle
!v4!              lla = llstr_dist_m_alph(jdist,mtot-mval_beta(istr))
!v4!              ula =  nstr_dist_m_alph(jdist,mtot-mval_beta(istr)) + lla - 1
!v4              !old
!v4              lla = lstr_alph_dist(1,jdist)
!v4              ula = lstr_alph_dist(2,jdist)
!v4              !old
!v4              do jstr = lla, ula
!v4                 !old
!v4                 if (mval_alph(jstr)+mval_beta(istr) /= mtot) cycle
!v4                 !old
!v4                 tmp = czero
!v4                 do j1x_m = 1, n1x_m_alph(m_kl,jstr)
!v4                    j1x = map1x_m_alph(j1x_m, m_kl, jstr)
!v4                    kfun = h1x_alph (j1x, jstr)
!v4                    lfun = p1x_alph (j1x, jstr)
!v4                    kstr = eq1x_alph(j1x, jstr)
!v4                    kdist = dist_str_alph(1,kstr)
!v4                    !old
!v4                    if (mval_alph(kstr)+mval_beta(lstr) /= mtot) cycle
!v4                    !old
!v4                    if (det_allowed(kdist,ldist) /= 0) then
!v4                       tmp = tmp + sgn1x_alph(j1x,jstr)*int2e(ifun,jfun,kfun,lfun)*cic(kstr,lstr)
!v4                    end if
!v4                 end do
!v4                 hcic(jstr,istr) = hcic(jstr,istr) + tmp * sgn1x_beta(i1x,istr)
!v4              end do
!v4           end do
!v4        end do
!v4     end do
!v4  end do
!v4  !$omp end do
!v4  !$omp end parallel

!v3  !$omp parallel default(shared) private(idist,m_kl,i1x,ifun,jfun,&
!v3  !$omp & lstr,ldist,lla,ula,j1x,kfun,lfun,kstr,kdist,tmp)
!v3  !$omp do
!v3  do istr = 1, nstr_beta
!v3     idist = dist_str_beta(1,istr)
!v3     do m_ij = -mmax2, mmax2
!v3        m_kl = -m_ij
!v3        do i1x_m = 1, n1x_m_beta(m_ij,istr)
!v3           i1x = map1x_m_beta(i1x_m,m_ij,istr)
!v3           ifun = h1x_beta  (i1x,istr)
!v3           jfun = p1x_beta  (i1x,istr)
!v3           lstr = eq1x_beta (i1x,istr)
!v3           ldist = dist_str_beta(1,lstr)
!v3           do jdist = 1, ndist_alph
!v3              if (det_allowed(jdist,idist) == 0) cycle
!v3              lla = llstr_dist_m_alph(jdist,mtot-mval_beta(istr))
!v3              ula =  nstr_dist_m_alph(jdist,mtot-mval_beta(istr)) + lla - 1
!v3              do jstr = lla, ula
!v3                 tmp = czero
!v3                 do j1x_m = 1, n1x_m_alph(m_kl,jstr)
!v3                    j1x = map1x_m_alph(j1x_m, m_kl, jstr)
!v3                    kfun = h1x_alph (j1x, jstr)
!v3                    lfun = p1x_alph (j1x, jstr)
!v3                    kstr = eq1x_alph(j1x, jstr)
!v3                    kdist = dist_str_alph(1,kstr)
!v3                    if (det_allowed(kdist,ldist) /= 0) then
!v3                       tmp = tmp + sgn1x_alph(j1x,jstr)*int2e(ifun,jfun,kfun,lfun)*cic(kstr,lstr)
!v3                    end if
!v3                 end do
!v3                 hcic(jstr,istr) = hcic(jstr,istr) + tmp * sgn1x_beta(i1x,istr)
!v3              end do
!v3           end do
!v3        end do
!v3     end do
!v3  end do
!v3  !$omp end do
!v3  !$omp end parallel

!v2  !$omp parallel default(shared) private(idist,m_kl,i1x,ifun,jfun,&
!v2  !$omp & lstr,ldist,lla,ula,j1x,kfun,lfun,kstr,kdist,tmp)
!v2  !$omp do
!v2  do istr = 1, nstr_beta
!v2     idist = dist_str_beta(1,istr)
!v2     do m_ij = -mmax2, mmax2
!v2        m_kl = -m_ij
!v2        do i1x_m = 1, n1x_m_beta(m_ij,istr)
!v2           i1x = map1x_m_beta(i1x_m,m_ij,istr)
!v2           ifun = h1x_beta  (i1x,istr)
!v2           jfun = p1x_beta  (i1x,istr)
!v2           lstr = eq1x_beta (i1x,istr)
!v2           ldist = dist_str_beta(1,lstr)
!v2           do jdist = 1, ndist_alph
!v2              lla = llstr_dist_m_alph(jdist,mtot-mval_beta(istr))
!v2              ula =  nstr_dist_m_alph(jdist,mtot-mval_beta(istr)) + lla - 1
!v2              do jstr = lla, ula
!v2                 tmp = czero
!v2                 do j1x_m = 1, n1x_m_alph(m_kl,jstr)
!v2                    j1x = map1x_m_alph(j1x_m, m_kl, jstr)
!v2                    kfun = h1x_alph (j1x, jstr)
!v2                    lfun = p1x_alph (j1x, jstr)
!v2                    kstr = eq1x_alph(j1x, jstr)
!v2                    kdist = dist_str_alph(1,kstr)
!v2                    if (det_allowed(kdist,ldist) /= 0) then
!v2                       tmp = tmp + sgn1x_alph(j1x,jstr)*int2e(ifun,jfun,kfun,lfun)*cic(kstr,lstr)
!v2                    end if
!v2                 end do
!v2                 hcic(jstr,istr) = hcic(jstr,istr) + tmp * sgn1x_beta(i1x,istr)
!v2              end do
!v2           end do
!v2        end do
!v2     end do
!v2  end do
!v2  !$omp end do
!v2  !$omp end parallel
end subroutine ormas_hcic_ras_abp
!######################################################################
!######################################################################
subroutine ormas_hcic_ras_effint1e(int1e, int2e, eff1e)

  use, intrinsic :: iso_c_binding
  use mod_const, only : czero, chalf
  use mod_ormas, only : nact

  implicit none
  complex(c_double_complex), intent(in) :: int1e(1:nact, 1:nact)
  complex(c_double_complex), intent(in) :: int2e(1:nact, 1:nact, 1:nact, 1:nact)
  complex(c_double_complex), intent(out) :: eff1e(1:nact, 1:nact)

  integer(c_long) :: iact, jact, kact
  complex(c_double_complex) :: tmp

  do iact = 1, nact
     do jact = 1, nact
        tmp = czero
        do kact = 1, jact - 1
           tmp = tmp + int2e(jact, kact, kact, iact)
        end do
        eff1e(jact, iact) = int1e(jact, iact) - tmp

        if (jact > iact) then
           eff1e(jact, iact) = eff1e(jact, iact) - int2e(jact, jact, jact, iact)
        end if
     end do
     eff1e(iact, iact) = eff1e(iact, iact) - int2e(iact, iact, iact, iact) * chalf
  end do

end subroutine ormas_hcic_ras_effint1e
!######################################################################
!######################################################################
subroutine ormas_hcic_ras_effint2e(int2e, eff2e)

  use, intrinsic :: iso_c_binding
  use mod_const, only : chalf
  use mod_ormas, only : nact

  implicit none
  complex(c_double_complex), intent(in) :: int2e(1:nact, 1:nact, 1:nact, 1:nact)
  complex(c_double_complex), intent(out) :: eff2e(1:nact, 1:nact, 1:nact, 1:nact)

  integer(c_long) :: iact, jact

  eff2e(1:nact, 1:nact, 1:nact, 1:nact) = int2e(1:nact, 1:nact, 1:nact, 1:nact)
  do iact = 1, nact
     do jact = 1, nact
        eff2e(jact, iact, jact, iact) = eff2e(jact, iact, jact, iact) * chalf
     end do
  end do

end subroutine ormas_hcic_ras_effint2e
!######################################################################
