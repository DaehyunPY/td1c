!######################################################################
subroutine ormas_hcic1_old(int1e, cic, hcic)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nact, nstr_alph, nstr_beta

  implicit none
  complex(c_double_complex), intent(in) :: int1e(1:nact, 1:nact)
  complex(c_double_complex), intent(in) :: cic(1:nstr_alph, 1:nstr_beta)
  complex(c_double_complex), intent(inout) :: hcic(1:nstr_alph, 1:nstr_beta)

  call ormas_hcic1_old_bbp(int1e, cic, hcic)
  call ormas_hcic1_old_aap(int1e, cic, hcic)

end subroutine ormas_hcic1_old
!######################################################################
subroutine ormas_hcic1_old_bbp(int1e, cic, hcic)

  use, intrinsic :: iso_c_binding
  use mod_sph, only : mmax2
  use mod_ormas, only : thrcic
  use mod_const, only : czero, chalf
  use mod_ormas, only : ncore, nact, nsub, lorb_sub, sub_orb, det_allowed
  use mod_ormas, only : ndist_alph, nstr_alph, lstr_alph_dist
  use mod_ormas, only : min_sub_beta, max_sub_beta, ndist_beta, dist_beta, &
       & nstr_beta, lstr_beta_dist, dist_str_beta, n1x_beta, p1x_beta, &
       & h1x_beta, eq1x_beta, sgn1x_beta
  use mod_ormas, only : n1x_m_alph, map1x_m_alph
  use mod_ormas, only : n1x_m_beta, map1x_m_beta

  implicit none
  complex(c_double_complex), intent(in) :: int1e(1:nact, 1:nact)
  complex(c_double_complex), intent(in) :: cic(1:nstr_alph, 1:nstr_beta)
  complex(c_double_complex), intent(inout) :: hcic(1:nstr_alph, 1:nstr_beta)

  complex(c_double_complex), allocatable :: sint(:,:)
  logical :: dotype2, dotype3
  integer(c_long) :: istr, jstr, kstr, lstr, i1x, j1x, ifun, jfun, kfun, lfun, &
       & idist, jdist, kdist, ii, lla, ula, llb, ulb, subi, subj, subk, subl, &
       & occi, occk, ksub, iproc, nproc, i1x_m, j1x_m, m_ij, m_kl
  complex(c_double_complex) :: tsgn, tmp, eff1e
  integer(c_long), external :: util_omp_nproc
  integer(c_long), external :: util_omp_iproc

  nproc = util_omp_nproc()
  allocate(sint(1:nstr_beta, 0:(nproc-1)))

!$omp parallel default(shared) private(idist,kstr,lstr,ifun,jfun,kfun,lfun, &
!$omp & lla,ula,llb,ulb,tsgn,tmp,subi,subj,eff1e,occi,occk,dotype2,dotype3, &
!$omp & subk,subl,i1x,j1x,i1x_m,j1x_m,m_ij,m_kl,iproc)
  iproc = util_omp_iproc()
!$omp do
  do istr = 1, nstr_beta
     idist = dist_str_beta(1, istr)
     call util_zcopy(nstr_beta, czero, 0, sint(1,iproc), 1)
!M-adapt
!    do i1x = 1, n1x_beta(0, istr)
     do m_ij = -mmax2, mmax2
        m_kl = -m_ij
        do i1x_m = 1, n1x_m_beta(m_ij, istr)
           i1x = map1x_m_beta(i1x_m, m_ij, istr)
!M-adapt
           ifun = h1x_beta  (i1x, istr)
           jfun = p1x_beta  (i1x, istr)
           kstr = eq1x_beta (i1x, istr)
           !TMP
           !if (kstr <= 0) cycle
           !TMP

           subi = sub_orb(ifun)
           subj = sub_orb(jfun)
           occi = dist_beta(subi, idist)
           tsgn = sgn1x_beta(i1x, istr) * chalf

           tmp = czero
!2e           if (occi > min_sub_beta(subi)) then
!2e              do ksub = 1, nsub
!2e                 if (ksub /= subi) then
!2e                    occk = dist_beta(ksub, idist)
!2e                    if (occk < max_sub_beta(ksub)) then
!2e                       do kfun = lorb_sub(1, ksub), lorb_sub(2, ksub)
!2e                          tmp = tmp + int2e(ifun, kfun, kfun, jfun)
!2e                       end do
!2e                    end if
!2e                 end if
!2e              end do
!2e           end if
           eff1e = int1e(ifun, jfun) - tmp * chalf
           sint(kstr, iproc) = sint(kstr, iproc) + sgn1x_beta(i1x, istr) * eff1e

!2e           dotype2 = (subi == subj) .and. (occi == min_sub_beta(subi))
!2e           dotype3 = (subi == subj) .and. (occi == max_sub_beta(subi))
!2e!NEW
!2e!          dotype2 = dotype2 .and. occi > 0
!2e!          dotype3 = dotype3 .and. occi < min(nelact(2), norb_sub(subi))
!2e!NEW
!2e
!2e!M-adapt
!2e!          do j1x = 1, n1x_beta(0, kstr)
!2e           do j1x_m = 1, n1x_m_beta(m_kl, kstr)
!2e              j1x = map1x_m_beta(j1x_m, m_kl, kstr)
!2e!M-adapt
!2e              kfun = h1x_beta (j1x, kstr)
!2e              lfun = p1x_beta (j1x, kstr)
!2e              lstr = eq1x_beta(j1x, kstr)
!2e              subk = sub_orb(kfun)
!2e              subl = sub_orb(lfun)
!2e              occk = dist_beta(subk, idist)
!2e
!2e!             if (mval(ncore+ifun) + mval(ncore+kfun) - mval(ncore+jfun) - mval(ncore+lfun) .ne. 0) cycle
!2e
!2e              tmp = int2e(ifun, jfun, kfun, lfun)
!2e              if (subi /= subk .and. subi /= subl) then
!2e                 if (dotype2) then
!2e                    tmp = tmp - int2e(ifun, lfun, kfun, jfun)
!2e                 end if
!2e                 if (dotype3 .and. occk > min_sub_beta(subk)) then
!2e                    tmp = tmp - int2e(kfun, jfun, ifun, lfun)
!2e                 end if
!2e              end if
!2e              sint(lstr, iproc) = sint(lstr, iproc) &
!2e                   & + tsgn * sgn1x_beta(j1x, kstr) * tmp
!2e           end do
        end do
     end do

     !WARNING: a room for improvement
     do jdist = 1, ndist_alph
        if (det_allowed(jdist, idist) == 0) cycle
        lla = lstr_alph_dist(1, jdist)
        ula = lstr_alph_dist(2, jdist)
        do kdist = 1, ndist_beta
           if (det_allowed(jdist, kdist) == 0) cycle
           llb = lstr_beta_dist(1, kdist)
           ulb = lstr_beta_dist(2, kdist)
           do jstr = llb, ulb
              if (abs(sint(jstr, iproc)) > thrcic) then
                 do ii = lla, ula
                    hcic(ii, istr) = hcic(ii, istr) + sint(jstr, iproc) * cic(ii, jstr)
                 end do
              end if
           end do
        end do
     end do
  end do
!$omp end do
!$omp end parallel

  deallocate(sint)

end subroutine ormas_hcic1_old_bbp
!######################################################################
!######################################################################
subroutine ormas_hcic1_old_aap(int1e, cic, hcic)

  use, intrinsic :: iso_c_binding
  use mod_sph, only : mmax2
  use mod_ormas, only : thrcic
  use mod_bas, only : smul
  use mod_const, only : czero, chalf
  use mod_ormas, only : ncore, nact, nsub, lorb_sub, sub_orb, det_allowed
  use mod_ormas, only : ndist_beta, nstr_beta, lstr_beta_dist
  use mod_ormas, only : min_sub_alph, max_sub_alph, ndist_alph, dist_alph, &
       & nstr_alph, lstr_alph_dist, dist_str_alph, n1x_alph, p1x_alph, &
       & h1x_alph, eq1x_alph, sgn1x_alph, nelact
  use mod_ormas, only : n1x_m_alph, map1x_m_alph
  use mod_ormas, only : n1x_m_beta, map1x_m_beta

  implicit none
  complex(c_double_complex), intent(in) :: int1e(1:nact, 1:nact)
  complex(c_double_complex), intent(in) :: cic(1:nstr_alph, 1:nstr_beta)
  complex(c_double_complex), intent(inout) :: hcic(1:nstr_alph, 1:nstr_beta)

  complex(c_double_complex), allocatable :: sint(:,:)
  logical :: dotype2, dotype3
  integer(c_long) :: istr, jstr, kstr, lstr, i1x, j1x, ifun, jfun, kfun, lfun, &
       & idist, jdist, kdist, ii, lla, ula, llb, ulb, subi, subj, subk, subl, &
       & occi, occk, ksub, iproc, nproc, i1x_m, j1x_m, m_ij, m_kl
  complex(c_double_complex) :: tsgn, tmp, eff1e
  integer(c_long), external :: util_omp_nproc
  integer(c_long), external :: util_omp_iproc

  nproc = util_omp_nproc()

  !##### restricted case #####
  if (smul == 1 .and. nelact(1) == nelact(2)) then
     do istr = 1, nstr_beta
        do jstr = 1, istr
           hcic(jstr, istr) = hcic(jstr, istr) + hcic(istr, jstr)
        end do
     end do

     do istr = 1, nstr_beta
        do jstr = istr + 1, nstr_beta
           hcic(jstr, istr) = hcic(istr, jstr)
        end do
     end do
     return
  end if


  !##### unrestricted case #####
  allocate(sint(1:nstr_alph, 0:(nproc-1)))

!$omp parallel default(shared) private(idist,kstr,lstr,ifun,jfun,kfun,lfun, &
!$omp & lla,ula,llb,ulb,tsgn,tmp,subi,subj,eff1e,occi,occk,dotype2,dotype3, &
!$omp & subk,subl,i1x,j1x,i1x_m,j1x_m,m_ij,m_kl,iproc)
  iproc = util_omp_iproc()
!$omp do
  do istr = 1, nstr_alph
     idist = dist_str_alph(1, istr)
     call util_zcopy(nstr_alph, czero, 0, sint(1,iproc), 1)
!M-adapt
!    do i1x = 1, n1x_alph(0, istr)
     do m_ij = -mmax2, mmax2
        m_kl = -m_ij
        do i1x_m = 1, n1x_m_alph(m_ij, istr)
           i1x = map1x_m_alph(i1x_m, m_ij, istr)
!M-adapt
           ifun = h1x_alph  (i1x, istr)
           jfun = p1x_alph  (i1x, istr)
           subi = sub_orb(ifun)
           subj = sub_orb(jfun)
           occi = dist_alph(subi, idist)
           kstr = eq1x_alph (i1x, istr)
           !TMP
           !if (kstr <= 0) cycle
           !TMP

           tsgn = sgn1x_alph(i1x, istr) * chalf

           tmp = czero
!2e           if (occi > min_sub_alph(subi)) then
!2e              do ksub = 1, nsub
!2e                 if (ksub /= subi) then
!2e                    occk = dist_alph(ksub, idist)
!2e                    if (occk < max_sub_alph(ksub)) then
!2e                       do kfun = lorb_sub(1, ksub), lorb_sub(2, ksub)
!2e                          tmp = tmp + int2e(ifun, kfun, kfun, jfun)
!2e                       end do
!2e                    end if
!2e                 end if
!2e              end do
!2e           end if
           eff1e = int1e(ifun, jfun) - tmp * chalf
           sint(kstr, iproc) = sint(kstr, iproc) &
                & + sgn1x_alph(i1x, istr) * eff1e

!2e           dotype2 = (subi == subj) .and. (occi == min_sub_alph(subi))
!2e           dotype3 = (subi == subj) .and. (occi == max_sub_alph(subi))
!2e!NEW
!2e!          dotype2 = dotype2 .and. occi > 0
!2e!          dotype3 = dotype3 .and. occi < min(nelact(2), norb_sub(subi))
!2e!NEW
!2e
!2e!M-adapt
!2e!          do j1x = 1, n1x_alph(0, kstr)
!2e           do j1x_m = 1, n1x_m_alph(m_kl, kstr)
!2e              j1x = map1x_m_alph(j1x_m, m_kl, kstr)
!2e!M-adapt
!2e              kfun = h1x_alph (j1x, kstr)
!2e              lfun = p1x_alph (j1x, kstr)
!2e              lstr = eq1x_alph(j1x, kstr)
!2e              subk = sub_orb(kfun)
!2e              subl = sub_orb(lfun)
!2e              occk = dist_alph(subk, idist)
!2e
!2e!             if (mval(ncore+ifun) + mval(ncore+kfun) - mval(ncore+jfun) - mval(ncore+lfun) .ne. 0) cycle
!2e
!2e              tmp = int2e(ifun, jfun, kfun, lfun)
!2e              if (subi /= subk .and. subi /= subl) then
!2e                 if (dotype2) then
!2e                    tmp = tmp - int2e(ifun, lfun, kfun, jfun)
!2e                 end if
!2e                 if (dotype3 .and. occk > min_sub_alph(subk)) then
!2e                    tmp = tmp - int2e(kfun, jfun, ifun, lfun)
!2e                 end if
!2e              end if
!2e              sint(lstr, iproc) = sint(lstr, iproc) &
!2e                   & + tsgn * sgn1x_alph(j1x, kstr) * tmp
!2e           end do
        end do
     end do

     !WARNING: a room for improvement
     do jdist = 1, ndist_beta
        if (det_allowed(idist, jdist) == 0) cycle
        llb = lstr_beta_dist(1, jdist)
        ulb = lstr_beta_dist(2, jdist)
        do kdist = 1, ndist_alph
           if (det_allowed(kdist, jdist) == 0) cycle
           lla = lstr_alph_dist(1, kdist)
           ula = lstr_alph_dist(2, kdist)
           do jstr = lla, ula
              if (abs(sint(jstr, iproc)) > thrcic) then
                 do ii = llb, ulb
                    hcic(istr, ii) = hcic(istr, ii) + sint(jstr, iproc) * cic(jstr, ii)
                 end do
              end if
           end do
        end do
     end do
  end do
!$omp end do
!$omp end parallel

  deallocate(sint)

end subroutine ormas_hcic1_old_aap
!######################################################################
