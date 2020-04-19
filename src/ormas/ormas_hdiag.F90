!######################################################################
subroutine ormas_hdiag(int1e, int2e, cic, hamd)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nact, nstr_alph, nstr_beta

  implicit none
  complex(c_double_complex), intent(in) :: int1e(1:nact, 1:nact)
  complex(c_double_complex), intent(in) :: int2e(1:nact, 1:nact, 1:nact, 1:nact)
  complex(c_double_complex), intent(in) :: cic(1:nstr_alph, 1:nstr_beta)
  complex(c_double_complex), intent(inout) :: hamd(1:nstr_alph, 1:nstr_beta)
  complex(c_double_complex), allocatable :: eff1e(:,:)
  stop "ormas_hdiag no longer supported."

!OLD  allocate(eff1e(1:nact, 1:nact))
!OLD
!OLD  ! 1e integrals with separable 2e contributions
!OLD  call ormas_hdiag_effint1e(int1e, int2e, eff1e)
!OLD
!OLD  ! NOTE: following order DOES matter
!OLD  call ormas_hdiag_bbp(eff1e, int2e, cic, hamd)
!OLD  call ormas_hdiag_aap(eff1e, int2e, cic, hamd)
!OLD  call ormas_hdiag_abp(int2e, cic, hamd)
!OLD
!OLD  deallocate(eff1e)

end subroutine ormas_hdiag
!######################################################################
!OLD!######################################################################
!OLDsubroutine ormas_hdiag_bbp(int1e, int2e, cic, hamd)
!OLD
!OLD  use, intrinsic :: iso_c_binding
!OLD  use mod_sph, only : mmax2
!OLD  use mod_ormas, only : thrcic
!OLD  use mod_const, only : czero, chalf
!OLD  use mod_ormas, only : ncore, nact, nsub, lorb_sub, sub_orb, det_allowed
!OLD  use mod_ormas, only : ndist_alph, nstr_alph, lstr_alph_dist
!OLD  use mod_ormas, only : min_sub_beta, max_sub_beta, ndist_beta, dist_beta, &
!OLD       & nstr_beta, lstr_beta_dist, dist_str_beta, n1x_beta, p1x_beta, &
!OLD       & h1x_beta, eq1x_beta, sgn1x_beta
!OLD  use mod_ormas, only : n1x_m_alph, map1x_m_alph
!OLD  use mod_ormas, only : n1x_m_beta, map1x_m_beta
!OLD
!OLD  implicit none
!OLD  complex(c_double_complex), intent(in) :: int1e(1:nact, 1:nact)
!OLD  complex(c_double_complex), intent(in) :: int2e(1:nact, 1:nact, 1:nact, 1:nact)
!OLD  complex(c_double_complex), intent(in) :: cic(1:nstr_alph, 1:nstr_beta)
!OLD  complex(c_double_complex), intent(inout) :: hamd(1:nstr_alph, 1:nstr_beta)
!OLD
!OLD  complex(c_double_complex), allocatable :: sint(:,:)
!OLD  logical :: dotype2, dotype3
!OLD  integer(c_int) :: istr, jstr, kstr, lstr, i1x, j1x, ifun, jfun, kfun, lfun, &
!OLD       & idist, jdist, kdist, ii, lla, ula, llb, ulb, subi, subj, subk, subl, &
!OLD       & occi, occk, ksub, iproc, nproc, i1x_m, j1x_m, m_ij, m_kl
!OLD  complex(c_double_complex) :: tsgn, tmp, eff1e
!OLD  integer(c_int), external :: util_omp_nproc
!OLD  integer(c_int), external :: util_omp_iproc
!OLD
!OLD  nproc = util_omp_nproc()
!OLD  allocate(sint(1:nstr_beta, 0:(nproc-1)))
!OLD
!OLD!$omp parallel default(shared) private(idist,kstr,lstr,ifun,jfun,kfun,lfun, &
!OLD!$omp & lla,ula,llb,ulb,tsgn,tmp,subi,subj,eff1e,occi,occk,dotype2,dotype3, &
!OLD!$omp & subk,subl,i1x,j1x,i1x_m,j1x_m,m_ij,m_kl,iproc)
!OLD  iproc = util_omp_iproc()
!OLD!$omp do
!OLD  do istr = 1, nstr_beta
!OLD     idist = dist_str_beta(1, istr)
!OLD     call util_zcopy(nstr_beta, czero, 0, sint(1,iproc), 1)
!OLD!M-adapt
!OLD!    do i1x = 1, n1x_beta(0, istr)
!OLD     do m_ij = -mmax2, mmax2
!OLD        m_kl = -m_ij
!OLD        do i1x_m = 1, n1x_m_beta(m_ij, istr)
!OLD           i1x = map1x_m_beta(i1x_m, m_ij, istr)
!OLD!M-adapt
!OLD           ifun = h1x_beta  (i1x, istr)
!OLD           jfun = p1x_beta  (i1x, istr)
!OLD           kstr = eq1x_beta (i1x, istr)
!OLD           subi = sub_orb(ifun)
!OLD           subj = sub_orb(jfun)
!OLD           occi = dist_beta(subi, idist)
!OLD           tsgn = sgn1x_beta(i1x, istr) * chalf
!OLD
!OLD           tmp = czero
!OLD           if (occi > min_sub_beta(subi)) then
!OLD              do ksub = 1, nsub
!OLD                 if (ksub /= subi) then
!OLD                    occk = dist_beta(ksub, idist)
!OLD                    if (occk < max_sub_beta(ksub)) then
!OLD                       do kfun = lorb_sub(1, ksub), lorb_sub(2, ksub)
!OLD                          tmp = tmp + int2e(ifun, kfun, kfun, jfun)
!OLD                       end do
!OLD                    end if
!OLD                 end if
!OLD              end do
!OLD           end if
!OLD           eff1e = int1e(ifun, jfun) - tmp * chalf
!OLD           sint(kstr, iproc) = sint(kstr, iproc) &
!OLD                & + sgn1x_beta(i1x, istr) * eff1e
!OLD
!OLD           dotype2 = (subi == subj) .and. (occi == min_sub_beta(subi))
!OLD           dotype3 = (subi == subj) .and. (occi == max_sub_beta(subi))
!OLD!NEW
!OLD!          dotype2 = dotype2 .and. occi > 0
!OLD!          dotype3 = dotype3 .and. occi < min(nelact(2), norb_sub(subi))
!OLD!NEW
!OLD
!OLD!M-adapt
!OLD!          do j1x = 1, n1x_beta(0, kstr)
!OLD           do j1x_m = 1, n1x_m_beta(m_kl, kstr)
!OLD              j1x = map1x_m_beta(j1x_m, m_kl, kstr)
!OLD!M-adapt
!OLD              kfun = h1x_beta (j1x, kstr)
!OLD              lfun = p1x_beta (j1x, kstr)
!OLD              lstr = eq1x_beta(j1x, kstr)
!OLD              subk = sub_orb(kfun)
!OLD              subl = sub_orb(lfun)
!OLD              occk = dist_beta(subk, idist)
!OLD
!OLD!             if (mval(ncore+ifun) + mval(ncore+kfun) - mval(ncore+jfun) - mval(ncore+lfun) .ne. 0) cycle
!OLD
!OLD              tmp = int2e(ifun, jfun, kfun, lfun)
!OLD              if (subi /= subk .and. subi /= subl) then
!OLD                 if (dotype2) then
!OLD                    tmp = tmp - int2e(ifun, lfun, kfun, jfun)
!OLD                 end if
!OLD                 if (dotype3 .and. occk > min_sub_beta(subk)) then
!OLD                    tmp = tmp - int2e(kfun, jfun, ifun, lfun)
!OLD                 end if
!OLD              end if
!OLD              sint(lstr, iproc) = sint(lstr, iproc) &
!OLD                   & + tsgn * sgn1x_beta(j1x, kstr) * tmp
!OLD           end do
!OLD        end do
!OLD     end do
!OLD
!OLD     !WARNING: a room for improvement
!OLD     do jdist = 1, ndist_alph
!OLD        if (det_allowed(jdist, idist) == 0) cycle
!OLD        lla = lstr_alph_dist(1, jdist)
!OLD        ula = lstr_alph_dist(2, jdist)
!OLD        do kdist = 1, ndist_beta
!OLD           if (det_allowed(jdist, kdist) == 0) cycle
!OLD           llb = lstr_beta_dist(1, kdist)
!OLD           ulb = lstr_beta_dist(2, kdist)
!OLD           do jstr = llb, ulb
!OLD              if (abs(sint(jstr, iproc)) > thrcic) then
!OLD                 do ii = lla, ula
!OLD                    hamd(ii, istr) = hamd(ii, istr) + sint(jstr, iproc) * cic(ii, jstr)
!OLD                 end do
!OLD              end if
!OLD           end do
!OLD        end do
!OLD     end do
!OLD  end do
!OLD!$omp end do
!OLD!$omp end parallel
!OLD
!OLD  deallocate(sint)
!OLD
!OLDend subroutine ormas_hdiag_bbp
!OLD!######################################################################
!OLD!######################################################################
!OLDsubroutine ormas_hdiag_aap(int1e, int2e, cic, hamd)
!OLD
!OLD  use, intrinsic :: iso_c_binding
!OLD  use mod_sph, only : mmax2
!OLD  use mod_ormas, only : thrcic
!OLD  use mod_bas, only : smul
!OLD  use mod_const, only : czero, chalf
!OLD  use mod_ormas, only : ncore, nact, nsub, lorb_sub, sub_orb, det_allowed
!OLD  use mod_ormas, only : ndist_beta, nstr_beta, lstr_beta_dist
!OLD  use mod_ormas, only : min_sub_alph, max_sub_alph, ndist_alph, dist_alph, &
!OLD       & nstr_alph, lstr_alph_dist, dist_str_alph, n1x_alph, p1x_alph, &
!OLD       & h1x_alph, eq1x_alph, sgn1x_alph, nelact
!OLD  use mod_ormas, only : n1x_m_alph, map1x_m_alph
!OLD  use mod_ormas, only : n1x_m_beta, map1x_m_beta
!OLD
!OLD  implicit none
!OLD  complex(c_double_complex), intent(in) :: int1e(1:nact, 1:nact)
!OLD  complex(c_double_complex), intent(in) :: int2e(1:nact, 1:nact, 1:nact, 1:nact)
!OLD  complex(c_double_complex), intent(in) :: cic(1:nstr_alph, 1:nstr_beta)
!OLD  complex(c_double_complex), intent(inout) :: hamd(1:nstr_alph, 1:nstr_beta)
!OLD
!OLD  complex(c_double_complex), allocatable :: sint(:,:)
!OLD  logical :: dotype2, dotype3
!OLD  integer(c_int) :: istr, jstr, kstr, lstr, i1x, j1x, ifun, jfun, kfun, lfun, &
!OLD       & idist, jdist, kdist, ii, lla, ula, llb, ulb, subi, subj, subk, subl, &
!OLD       & occi, occk, ksub, iproc, nproc, i1x_m, j1x_m, m_ij, m_kl
!OLD  complex(c_double_complex) :: tsgn, tmp, eff1e
!OLD  integer(c_int), external :: util_omp_nproc
!OLD  integer(c_int), external :: util_omp_iproc
!OLD
!OLD  nproc = util_omp_nproc()
!OLD
!OLD  !##### restricted case #####
!OLD  if (smul == 1 .and. nelact(1) == nelact(2)) then
!OLD     do istr = 1, nstr_beta
!OLD        do jstr = 1, istr
!OLD           hamd(jstr, istr) = hamd(jstr, istr) + hamd(istr, jstr)
!OLD        end do
!OLD     end do
!OLD
!OLD     do istr = 1, nstr_beta
!OLD        do jstr = istr + 1, nstr_beta
!OLD           hamd(jstr, istr) = hamd(istr, jstr)
!OLD        end do
!OLD     end do
!OLD     return
!OLD  end if
!OLD
!OLD
!OLD  !##### unrestricted case #####
!OLD  allocate(sint(1:nstr_alph, 0:(nproc-1)))
!OLD
!OLD!$omp parallel default(shared) private(idist,kstr,lstr,ifun,jfun,kfun,lfun, &
!OLD!$omp & lla,ula,llb,ulb,tsgn,tmp,subi,subj,eff1e,occi,occk,dotype2,dotype3, &
!OLD!$omp & subk,subl,i1x,j1x,i1x_m,j1x_m,m_ij,m_kl,iproc)
!OLD  iproc = util_omp_iproc()
!OLD!$omp do
!OLD  do istr = 1, nstr_alph
!OLD     idist = dist_str_alph(1, istr)
!OLD     call util_zcopy(nstr_alph, czero, 0, sint(1,iproc), 1)
!OLD!M-adapt
!OLD!    do i1x = 1, n1x_alph(0, istr)
!OLD     do m_ij = -mmax2, mmax2
!OLD        m_kl = -m_ij
!OLD        do i1x_m = 1, n1x_m_alph(m_ij, istr)
!OLD           i1x = map1x_m_alph(i1x_m, m_ij, istr)
!OLD!M-adapt
!OLD           ifun = h1x_alph  (i1x, istr)
!OLD           jfun = p1x_alph  (i1x, istr)
!OLD           subi = sub_orb(ifun)
!OLD           subj = sub_orb(jfun)
!OLD           occi = dist_alph(subi, idist)
!OLD           kstr = eq1x_alph (i1x, istr)
!OLD           tsgn = sgn1x_alph(i1x, istr) * chalf
!OLD
!OLD           tmp = czero
!OLD           if (occi > min_sub_alph(subi)) then
!OLD              do ksub = 1, nsub
!OLD                 if (ksub /= subi) then
!OLD                    occk = dist_alph(ksub, idist)
!OLD                    if (occk < max_sub_alph(ksub)) then
!OLD                       do kfun = lorb_sub(1, ksub), lorb_sub(2, ksub)
!OLD                          tmp = tmp + int2e(ifun, kfun, kfun, jfun)
!OLD                       end do
!OLD                    end if
!OLD                 end if
!OLD              end do
!OLD           end if
!OLD           eff1e = int1e(ifun, jfun) - tmp * chalf
!OLD           sint(kstr, iproc) = sint(kstr, iproc) &
!OLD                & + sgn1x_alph(i1x, istr) * eff1e
!OLD
!OLD           dotype2 = (subi == subj) .and. (occi == min_sub_alph(subi))
!OLD           dotype3 = (subi == subj) .and. (occi == max_sub_alph(subi))
!OLD!NEW
!OLD!          dotype2 = dotype2 .and. occi > 0
!OLD!          dotype3 = dotype3 .and. occi < min(nelact(2), norb_sub(subi))
!OLD!NEW
!OLD
!OLD!M-adapt
!OLD!          do j1x = 1, n1x_alph(0, kstr)
!OLD           do j1x_m = 1, n1x_m_alph(m_kl, kstr)
!OLD              j1x = map1x_m_alph(j1x_m, m_kl, kstr)
!OLD!M-adapt
!OLD              kfun = h1x_alph (j1x, kstr)
!OLD              lfun = p1x_alph (j1x, kstr)
!OLD              lstr = eq1x_alph(j1x, kstr)
!OLD              subk = sub_orb(kfun)
!OLD              subl = sub_orb(lfun)
!OLD              occk = dist_alph(subk, idist)
!OLD
!OLD!             if (mval(ncore+ifun) + mval(ncore+kfun) - mval(ncore+jfun) - mval(ncore+lfun) .ne. 0) cycle
!OLD
!OLD              tmp = int2e(ifun, jfun, kfun, lfun)
!OLD              if (subi /= subk .and. subi /= subl) then
!OLD                 if (dotype2) then
!OLD                    tmp = tmp - int2e(ifun, lfun, kfun, jfun)
!OLD                 end if
!OLD                 if (dotype3 .and. occk > min_sub_alph(subk)) then
!OLD                    tmp = tmp - int2e(kfun, jfun, ifun, lfun)
!OLD                 end if
!OLD              end if
!OLD              sint(lstr, iproc) = sint(lstr, iproc) &
!OLD                   & + tsgn * sgn1x_alph(j1x, kstr) * tmp
!OLD           end do
!OLD        end do
!OLD     end do
!OLD
!OLD     !WARNING: a room for improvement
!OLD     do jdist = 1, ndist_beta
!OLD        if (det_allowed(idist, jdist) == 0) cycle
!OLD        llb = lstr_beta_dist(1, jdist)
!OLD        ulb = lstr_beta_dist(2, jdist)
!OLD        do kdist = 1, ndist_alph
!OLD           if (det_allowed(kdist, jdist) == 0) cycle
!OLD           lla = lstr_alph_dist(1, kdist)
!OLD           ula = lstr_alph_dist(2, kdist)
!OLD           do jstr = lla, ula
!OLD              if (abs(sint(jstr, iproc)) > thrcic) then
!OLD                 do ii = llb, ulb
!OLD                    hamd(istr, ii) = hamd(istr, ii) + sint(jstr, iproc) * cic(jstr, ii)
!OLD                 end do
!OLD              end if
!OLD           end do
!OLD        end do
!OLD     end do
!OLD  end do
!OLD!$omp end do
!OLD!$omp end parallel
!OLD
!OLD  deallocate(sint)
!OLD
!OLDend subroutine ormas_hdiag_aap
!OLD!######################################################################
!OLD!######################################################################
!OLDsubroutine ormas_hdiag_abp(int2e, cic, hamd)
!OLD
!OLD  use, intrinsic :: iso_c_binding
!OLD  use mod_bas, only : mval
!OLD  use mod_sph, only : mmax2
!OLD  use mod_ormas, only : thrcic
!OLD  use mod_const, only : czero
!OLD  use mod_ormas, only : ncore, nact, det_allowed
!OLD  use mod_ormas, only : nstr_alph, &
!OLD       & dist_str_alph, n1xr_alph, l1xr_alph, r1xr_alph, sgn1xr_alph
!OLD  use mod_ormas, only : nstr_beta, &
!OLD       & dist_str_beta, n1x_beta, p1x_beta, h1x_beta, eq1x_beta, sgn1x_beta
!OLD  use mod_ormas, only : n1x_m_alph, map1x_m_alph
!OLD  use mod_ormas, only : n1x_m_beta, map1x_m_beta
!OLD
!OLD  implicit none
!OLD  complex(c_double_complex), intent(in) :: int2e(1:nact, 1:nact, 1:nact, 1:nact)
!OLD  complex(c_double_complex), intent(in) :: cic(1:nstr_alph, 1:nstr_beta)
!OLD  complex(c_double_complex), intent(inout) :: hamd(1:nstr_alph, 1:nstr_beta)
!OLD
!OLD  complex(c_double_complex), allocatable :: sint(:,:)
!OLD  complex(c_double_complex), allocatable :: scic(:,:)
!OLD  complex(c_double_complex), allocatable :: hscic(:,:)
!OLD  integer(c_int) :: istr, jstr, kstr, j1x, ifun, jfun, kfun, lfun, &
!OLD       & idist, jdist, kdist, ii, n1xra, iproc, nproc, j1x_m, m_ij, m_kl
!OLD  integer(c_int), external :: util_omp_nproc
!OLD  integer(c_int), external :: util_omp_iproc
!OLD
!OLD  nproc = util_omp_nproc()
!OLD  allocate(scic(1:nstr_alph, 1:nstr_beta))
!OLD  allocate(sint(1:nstr_beta, 0:(nproc-1)))
!OLD  allocate(hscic(1:nstr_alph, 0:(nproc-1)))
!OLD
!OLD  do ifun = 1, nact
!OLD     do jfun = 1, nact
!OLD
!OLD!$omp parallel default(shared) private(istr,idist,kdist,n1xra,iproc)
!OLD        iproc = util_omp_iproc()
!OLD!$omp do
!OLD        do kstr = 1, nstr_beta
!OLD           kdist = dist_str_beta(1, kstr)
!OLD           n1xra = n1xr_alph(ifun, jfun)
!OLD           call util_zcopy(n1xra, czero, 0, scic(1, kstr), 1)
!OLD           do ii = 1, n1xra
!OLD              istr = l1xr_alph(ii, ifun, jfun)
!OLD              idist = dist_str_alph(1, istr)
!OLD              if (det_allowed(idist, kdist) /= 0) then
!OLD                 scic(ii, kstr) = cic(istr, kstr) * sgn1xr_alph(ii, ifun, jfun)
!OLD              end if
!OLD           end do
!OLD        end do
!OLD!$omp end do
!OLD!$omp end parallel
!OLD
!OLD!$omp parallel default(shared) private(idist,jdist,istr,kstr,kfun,lfun,n1xra,j1x,j1x_m,m_ij,m_kl,iproc)
!OLD        iproc = util_omp_iproc()
!OLD        m_ij = mval(ncore+ifun) - mval(ncore+jfun)
!OLD        m_kl = -m_ij
!OLD!$omp do
!OLD        do jstr = 1, nstr_beta
!OLD           jdist = dist_str_beta(1, jstr)
!OLD           call util_zcopy(nstr_beta, czero, 0, sint(1,iproc), 1)
!OLD!M-adapt
!OLD!          do j1x = 1, n1x_beta(0, jstr)
!OLD           do j1x_m = 1, n1x_m_beta(m_kl, jstr)
!OLD              j1x = map1x_m_beta(j1x_m, m_kl, jstr)
!OLD!M-adapt
!OLD              kfun = h1x_beta (j1x, jstr)
!OLD              lfun = p1x_beta (j1x, jstr)
!OLD              kstr = eq1x_beta(j1x, jstr)
!OLD              sint(kstr, iproc) = sint(kstr, iproc) &
!OLD                   & + sgn1x_beta(j1x, jstr) * int2e(ifun, jfun, kfun, lfun)
!OLD           end do
!OLD
!OLD           n1xra = n1xr_alph(ifun, jfun)           
!OLD           call util_zcopy(n1xra, czero, 0, hscic(1, iproc), 1)
!OLD           do kstr = 1, nstr_beta
!OLD              if (abs(sint(kstr, iproc)) > thrcic) then
!OLD!                call zaxpy(n1xra, sint(kstr,iproc), scic(1, kstr), 1, hscic(1, iproc), 1)
!OLD                 do ii = 1, n1xra
!OLD                    hscic(ii, iproc) = hscic(ii, iproc) &
!OLD                         & + sint(kstr, iproc) * scic(ii, kstr)
!OLD                 end do
!OLD              end if
!OLD           end do
!OLD
!OLD           ! vectolized scattering
!OLD           do ii = 1, n1xra
!OLD              istr = r1xr_alph(ii, ifun, jfun)
!OLD              idist = dist_str_alph(1, istr)
!OLD              if (det_allowed(idist, jdist) /= 0) then
!OLD                 hamd(istr, jstr) = hamd(istr, jstr) + hscic(ii, iproc)
!OLD              end if
!OLD           end do
!OLD        end do
!OLD!$omp end do
!OLD!$omp end parallel
!OLD     end do
!OLD  end do
!OLD
!OLD  deallocate(hscic)
!OLD  deallocate(sint)
!OLD  deallocate(scic)
!OLD
!OLDend subroutine ormas_hdiag_abp
!OLD!######################################################################
!OLDsubroutine ormas_hdiag_effint1e(int1e, int2e, eff1e)
!OLD
!OLD  use, intrinsic :: iso_c_binding
!OLD  use mod_const, only : czero, chalf
!OLD  use mod_ormas, only : nact, nsub, lorb_sub
!OLD
!OLD  implicit none
!OLD  complex(c_double_complex), intent(in) :: int1e(1:nact, 1:nact)
!OLD  complex(c_double_complex), intent(in) :: int2e(1:nact, 1:nact, 1:nact, 1:nact)
!OLD  complex(c_double_complex), intent(out) :: eff1e(1:nact, 1:nact)
!OLD
!OLD  integer(c_int) :: iact, jact, kact, jsub
!OLD  complex(c_double_complex) :: tmp
!OLD
!OLD!debug  write(6, "('WARNING: skip 2e part in effint1e!')")
!OLD  do iact = 1, nact
!OLD     do jsub = 1, nsub
!OLD        do jact = lorb_sub(1, jsub), lorb_sub(2, jsub)
!OLD           tmp = czero
!OLD           do kact = lorb_sub(1, jsub), lorb_sub(2, jsub)
!OLD              tmp = tmp + int2e(jact, kact, kact, iact)
!OLD           end do
!OLD           eff1e(jact, iact) = int1e(jact, iact) - chalf * tmp
!OLD        end do
!OLD     end do
!OLD  end do
!OLD
!OLDend subroutine ormas_hdiag_effint1e
!OLD!######################################################################
