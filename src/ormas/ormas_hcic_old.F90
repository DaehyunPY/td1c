!######################################################################
subroutine ormas_hcic_old(int1e, int2e, cic, hcic)

  use, intrinsic :: iso_c_binding
  use mod_control, only : name
  use mod_ormas, only : nact, nstr_alph, nstr_beta

  implicit none
  complex(c_double_complex), intent(in) :: int1e(1:nact, 1:nact)
  complex(c_double_complex), intent(in) :: int2e(1:nact, 1:nact, 1:nact, 1:nact)
  complex(c_double_complex), intent(in) :: cic(1:nstr_alph, 1:nstr_beta)
  complex(c_double_complex), intent(inout) :: hcic(1:nstr_alph, 1:nstr_beta)
  complex(c_double_complex), allocatable :: eff1e(:,:), eff2e(:,:,:,:)

  stop "ormas_hcic_old no longer supported."

!OLD  allocate(eff1e(1:nact, 1:nact))
!OLD  allocate(eff2e(1:nact, 1:nact, 1:nact, 1:nact))
!OLD
!OLD  ! 1e integrals with separable 2e contributions
!OLD  call ormas_hcic_old_effint1e(int1e, int2e, eff1e)
!OLD  call ormas_hcic_old_effint2e(int2e, eff2e)
!OLD
!OLD  call ormas_hcic_old_bbp(eff1e, eff2e, cic, hcic)
!OLD  call ormas_hcic_old_aap(eff1e, eff2e, cic, hcic)
!OLD  call ormas_hcic_old_abp(int2e, cic, hcic)
!OLD
!OLD!debug
!OLD!  call ormas_cic_print(hcic, trim(name)//".hcic")
!OLD!  stop "STOP for debug @ ormas_hcic_old"
!OLD!debug
!OLD  deallocate(eff2e)
!OLD  deallocate(eff1e)

end subroutine ormas_hcic_old
!OLD!######################################################################
!OLDsubroutine ormas_hcic_old_bbp(int1e, int2e, cic, hcic)
!OLD
!OLD  use, intrinsic :: iso_c_binding
!OLD  use mod_ormas, only : thrcic
!OLD  use mod_const, only : czero
!OLD  use mod_bas, only : mtot
!OLD  use mod_sph, only : mmax2
!OLD  use mod_ormas, only : nact, det_allowed
!OLD  use mod_ormas, only : ndist_alph, nstr_alph, lstr_alph_dist
!OLD  use mod_ormas, only : ndist_beta, nstr_beta, lstr_beta_dist, &
!OLD       & dist_str_beta, n1x_beta, p1x_beta, h1x_beta, eq1x_beta, sgn1x_beta
!OLD  use mod_ormas, only : n1x_m_alph, map1x_m_alph
!OLD  use mod_ormas, only : n1x_m_beta, map1x_m_beta
!OLD  use mod_ormas, only : nstr_dist_m_alph, llstr_dist_m_alph, mval_alph
!OLD  use mod_ormas, only : nstr_dist_m_beta, llstr_dist_m_beta, mval_beta
!OLD
!OLD  implicit none
!OLD  complex(c_double_complex), intent(in) :: int1e(1:nact, 1:nact)
!OLD  complex(c_double_complex), intent(in) :: int2e(1:nact, 1:nact, 1:nact, 1:nact)
!OLD  complex(c_double_complex), intent(in) :: cic(1:nstr_alph, 1:nstr_beta)
!OLD  complex(c_double_complex), intent(inout) :: hcic(1:nstr_alph, 1:nstr_beta)
!OLD
!OLD  complex(c_double_complex), allocatable :: sint(:,:)
!OLD  integer(c_int) :: istr, jstr, kstr, lstr, i1x, j1x, ifun, jfun, kfun, lfun, &
!OLD       & iord, jord, idist, jdist, kdist, ii, lla, ula, llb, ulb, iproc, nproc, &
!OLD       & i1x_m, j1x_m, m_ij, m_kl
!OLD  complex(c_double_complex) :: tsgn
!OLD  integer(c_int), external :: util_omp_nproc
!OLD  integer(c_int), external :: util_omp_iproc
!OLD
!OLD  nproc = util_omp_nproc()
!OLD  allocate(sint(1:nstr_beta, 0:(nproc-1)))
!OLD
!OLD  !$omp parallel default(shared) &
!OLD  !$omp private(iproc,idist,kstr,lstr,ifun,jfun,kfun,lfun) &
!OLD  !$omp private(iord,jord,lla,ula,llb,ulb,tsgn,i1x,j1x,m_ij,m_kl)
!OLD  iproc = util_omp_iproc()
!OLD!$omp do
!OLD  do istr = 1, nstr_beta
!OLD     idist = dist_str_beta(1, istr)
!OLD     call util_zcopy(nstr_beta, czero, 0, sint(1,iproc), 1)
!OLD!M   do i1x = 1, n1x_beta(0, istr)
!OLD     do m_ij = -mmax2, mmax2
!OLD        m_kl = -m_ij
!OLD        do i1x_m = 1, n1x_m_beta(m_ij, istr)
!OLD           i1x = map1x_m_beta(i1x_m, m_ij, istr)
!OLD           ifun = h1x_beta  (i1x, istr)
!OLD           jfun = p1x_beta  (i1x, istr)
!OLD           kstr = eq1x_beta (i1x, istr)
!OLD           tsgn = sgn1x_beta(i1x, istr)
!OLD           iord = nact * (ifun - 1) + jfun
!OLD           if (m_ij == 0) sint(kstr, iproc) = sint(kstr, iproc) + sgn1x_beta(i1x, istr) * int1e(ifun, jfun)
!OLD!M         do j1x = 1, n1x_beta(0, kstr)
!OLD           do j1x_m = 1, n1x_m_beta(m_kl, kstr)
!OLD              j1x = map1x_m_beta(j1x_m, m_kl, kstr)
!OLD              kfun = h1x_beta (j1x, kstr)
!OLD              lfun = p1x_beta (j1x, kstr)
!OLD              lstr = eq1x_beta(j1x, kstr)
!OLD              jord = nact * (kfun - 1) + lfun
!OLD              if (jord <= iord) then
!OLD                 sint(lstr, iproc) = sint(lstr, iproc) + tsgn * sgn1x_beta(j1x, kstr) * int2e(ifun, jfun, kfun, lfun)
!OLD              end if
!OLD           end do
!OLD        end do
!OLD     end do
!OLD
!OLD!debug     do kdist = 1, ndist_beta
!OLD!debug        llb = lstr_beta_dist(1, kdist)
!OLD!debug        ulb = lstr_beta_dist(2, kdist)
!OLD!debug        do jdist = 1, ndist_alph
!OLD!debug!debug           write(6, "('<JI|JK>:', 4i4, l4)") jdist-1, idist-1, jdist-1, kdist-1, &
!OLD!debug!debug                (det_allowed(jdist,kdist)/=0).and.(det_allowed(jdist,idist)/=0)
!OLD!debug           if (det_allowed(jdist, kdist) == 0 .or. &
!OLD!debug               det_allowed(jdist, idist) == 0) cycle
!OLD!debug           lla = lstr_alph_dist(1, jdist)
!OLD!debug           ula = lstr_alph_dist(2, jdist)
!OLD!debug           do jstr = llb, ulb
!OLD!debug              if (abs(sint(jstr, iproc)) > thrcic) then
!OLD!debug                 do ii = lla, ula
!OLD!debug                    hcic(ii, istr) = hcic(ii, istr) + sint(jstr, iproc) * cic(ii, jstr)
!OLD!debug                 end do
!OLD!debug              end if
!OLD!debug           end do
!OLD!debug        end do
!OLD!debug     end do
!OLD
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
!OLD                    hcic(ii, istr) = hcic(ii, istr) + sint(jstr, iproc) * cic(ii, jstr)
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
!OLDend subroutine ormas_hcic_old_bbp
!OLD!######################################################################
!OLD!######################################################################
!OLDsubroutine ormas_hcic_old_aap(int1e, int2e, cic, hcic)
!OLD
!OLD  use, intrinsic :: iso_c_binding
!OLD  use mod_bas, only : smul
!OLD  use mod_ormas, only : thrcic
!OLD  use mod_const, only : czero
!OLD  use mod_ormas, only : nact, det_allowed, nelact
!OLD  use mod_ormas, only : ndist_alph, nstr_alph, lstr_alph_dist, &
!OLD       & dist_str_alph, n1x_alph, p1x_alph, h1x_alph, eq1x_alph, sgn1x_alph
!OLD  use mod_ormas, only : ndist_beta, nstr_beta, lstr_beta_dist
!OLD
!OLD  implicit none
!OLD  complex(c_double_complex), intent(in) :: int1e(1:nact, 1:nact)
!OLD  complex(c_double_complex), intent(in) :: int2e(1:nact, 1:nact, 1:nact, 1:nact)
!OLD  complex(c_double_complex), intent(in) :: cic(1:nstr_alph, 1:nstr_beta)
!OLD  complex(c_double_complex), intent(inout) :: hcic(1:nstr_alph, 1:nstr_beta)
!OLD
!OLD  complex(c_double_complex), allocatable :: sint(:,:)
!OLD  integer(c_int) :: istr, jstr, kstr, lstr, i1x, j1x, ifun, jfun, kfun, lfun, &
!OLD       & iord, jord, idist, jdist, kdist, ii, lla, ula, llb, ulb, iproc, nproc
!OLD  complex(c_double_complex) :: tsgn
!OLD  integer(c_int), external :: util_omp_nproc
!OLD  integer(c_int), external :: util_omp_iproc
!OLD
!OLD  nproc = util_omp_nproc()
!OLD
!OLD  !##### restricted case #####
!OLD  if (smul == 1 .and. nelact(1) == nelact(2)) then
!OLD     do istr = 1, nstr_beta
!OLD        do jstr = 1, istr
!OLD           hcic(jstr, istr) = hcic(jstr, istr) + hcic(istr, jstr)
!OLD        end do
!OLD     end do
!OLD     do istr = 1, nstr_beta
!OLD        do jstr = istr + 1, nstr_alph
!OLD           hcic(jstr, istr) = hcic(istr, jstr)
!OLD        end do
!OLD     end do
!OLD     return
!OLD  end if
!OLD
!OLD  !##### unrestricted case #####
!OLD  allocate(sint(1:nstr_alph, 0:(nproc-1)))
!OLD
!OLD!$omp parallel default(shared) private(iproc,idist,kstr,lstr,ifun,jfun,kfun,lfun,iord,jord,lla,ula,llb,ulb,tsgn)
!OLD  iproc = util_omp_iproc()
!OLD!$omp do
!OLD  do istr = 1, nstr_alph
!OLD     idist = dist_str_alph(1, istr)
!OLD     call util_zcopy(nstr_alph, czero, 0, sint(1,iproc), 1)
!OLD     do i1x = 1, n1x_alph(0, istr)
!OLD        ifun = h1x_alph  (i1x, istr)
!OLD        jfun = p1x_alph  (i1x, istr)
!OLD        kstr = eq1x_alph (i1x, istr)
!OLD        tsgn = sgn1x_alph(i1x, istr)
!OLD        iord = nact * (ifun - 1) + jfun
!OLD        sint(kstr, iproc) = sint(kstr, iproc) &
!OLD             & + sgn1x_alph(i1x, istr) * int1e(ifun, jfun)
!OLD
!OLD        do j1x = 1, n1x_alph(0, kstr)
!OLD           kfun = h1x_alph (j1x, kstr)
!OLD           lfun = p1x_alph (j1x, kstr)
!OLD           lstr = eq1x_alph(j1x, kstr)
!OLD           jord = nact * (kfun - 1) + lfun
!OLD           if (jord <= iord) then
!OLD              sint(lstr, iproc) = sint(lstr, iproc) &
!OLD                   & + tsgn * sgn1x_alph(j1x, kstr) * int2e(ifun, jfun, kfun, lfun)
!OLD           end if
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
!OLD                    hcic(istr, ii) = hcic(istr, ii) + sint(jstr, iproc) * cic(jstr, ii)
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
!OLDend subroutine ormas_hcic_old_aap
!OLD!######################################################################
!OLD!######################################################################
!OLD!######################################################################
!OLDsubroutine ormas_hcic_old_abp(int2e, cic, hcic)
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
!OLD  complex(c_double_complex), intent(inout) :: hcic(1:nstr_alph, 1:nstr_beta)
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
!OLD              !TMP
!OLD              !if (kstr <= 0) cycle
!OLD              !TMP
!OLD
!OLD              sint(kstr, iproc) = sint(kstr, iproc) &
!OLD                   + sgn1x_beta(j1x, jstr) * int2e(ifun, jfun, kfun, lfun)
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
!OLD                 hcic(istr, jstr) = hcic(istr, jstr) + hscic(ii, iproc)
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
!OLDend subroutine ormas_hcic_old_abp
!OLD!######################################################################
!OLD!######################################################################
!OLDsubroutine ormas_hcic_old_effint1e(int1e, int2e, eff1e)
!OLD
!OLD  use, intrinsic :: iso_c_binding
!OLD  use mod_const, only : czero, chalf
!OLD  use mod_ormas, only : nact
!OLD
!OLD  implicit none
!OLD  complex(c_double_complex), intent(in) :: int1e(1:nact, 1:nact)
!OLD  complex(c_double_complex), intent(in) :: int2e(1:nact, 1:nact, 1:nact, 1:nact)
!OLD  complex(c_double_complex), intent(out) :: eff1e(1:nact, 1:nact)
!OLD
!OLD  integer(c_int) :: iact, jact, kact
!OLD  complex(c_double_complex) :: tmp
!OLD
!OLD  do iact = 1, nact
!OLD     do jact = 1, nact
!OLD        tmp = czero
!OLD        do kact = 1, jact - 1
!OLD           tmp = tmp + int2e(jact, kact, kact, iact)
!OLD        end do
!OLD        eff1e(jact, iact) = int1e(jact, iact) - tmp
!OLD
!OLD        if (jact > iact) then
!OLD           eff1e(jact, iact) = eff1e(jact, iact) - int2e(jact, jact, jact, iact)
!OLD        end if
!OLD     end do
!OLD     eff1e(iact, iact) = eff1e(iact, iact) - int2e(iact, iact, iact, iact) * chalf
!OLD  end do
!OLD
!OLDend subroutine ormas_hcic_old_effint1e
!OLD!######################################################################
!OLD!######################################################################
!OLDsubroutine ormas_hcic_old_effint2e(int2e, eff2e)
!OLD
!OLD  use, intrinsic :: iso_c_binding
!OLD  use mod_const, only : chalf
!OLD  use mod_ormas, only : nact
!OLD
!OLD  implicit none
!OLD  complex(c_double_complex), intent(in) :: int2e(1:nact, 1:nact, 1:nact, 1:nact)
!OLD  complex(c_double_complex), intent(out) :: eff2e(1:nact, 1:nact, 1:nact, 1:nact)
!OLD
!OLD  integer(c_int) :: iact, jact
!OLD
!OLD  eff2e(1:nact, 1:nact, 1:nact, 1:nact) = int2e(1:nact, 1:nact, 1:nact, 1:nact)
!OLD  do iact = 1, nact
!OLD     do jact = 1, nact
!OLD        eff2e(jact, iact, jact, iact) = eff2e(jact, iact, jact, iact) * chalf
!OLD     end do
!OLD  end do
!OLD
!OLDend subroutine ormas_hcic_old_effint2e
!OLD!######################################################################
