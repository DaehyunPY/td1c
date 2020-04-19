!######################################################################
subroutine ormas_mkden2x_old(cic, den2x)

  use, intrinsic :: iso_c_binding
  use mod_const, only : czero
  use mod_ormas, only : s2zero, nact, nstr_alph, nstr_beta, nelact

  implicit none
  complex(c_double_complex), intent(in) :: cic(1:nstr_alph, 1:nstr_beta)
  complex(c_double_complex), intent(out) :: den2x(1:nact, 1:nact, 1:nact, 1:nact)
  complex(c_double_complex), allocatable :: den2xp(:,:,:,:,:)
  integer(c_long) :: nproc,iproc,iact,jact,kact,lact
  integer(c_long), external :: util_omp_nproc

  nproc = util_omp_nproc()
  allocate(den2xp(1:nact, 1:nact, 1:nact, 1:nact, 0:(nproc-1)))
  den2x(1:nact, 1:nact, 1:nact, 1:nact) = czero
  den2xp(1:nact, 1:nact, 1:nact, 1:nact, 0:(nproc-1)) = czero

  call ormas_mkden2x_old_bbp(nproc, cic, den2xp)
  call ormas_mkden2x_old_abp(nproc, cic, den2xp)
  if (nelact(1) == nelact(2)) then
     den2xp(1:nact, 1:nact, 1:nact, 1:nact, 0:(nproc-1)) = &
     den2xp(1:nact, 1:nact, 1:nact, 1:nact, 0:(nproc-1)) + &
     den2xp(1:nact, 1:nact, 1:nact, 1:nact, 0:(nproc-1))
  else
     stop 'ormas_mkden2x: only for S=0 and Ms=0'
!    call ormas_mkden2x_old_aap(cic, lwork, work, den2xp)
!    call ormas_mkden2x_old_bap(cic, lwork, work, den2xp)
  end if

  do iproc = 0, nproc - 1
     den2x (1:nact, 1:nact, 1:nact, 1:nact) = & ! WATCH the Sign! See BELOW
     den2x (1:nact, 1:nact, 1:nact, 1:nact) + & ! WATCH the Sign! See BELOW
     den2xp(1:nact, 1:nact, 1:nact, 1:nact, iproc)
  end do
  deallocate(den2xp)

!debug
!  write(6, "('ormas_mkden2x_old: den2x')")
!  do lact = 1, nact
!  do kact = 1, nact
!  do jact = 1, nact
!  do iact = 1, nact
!     write(6, "(2f20.10)") den2x(iact,jact,kact,lact)
!  end do
!  end do
!  end do
!  end do
!  stop
!debug

end subroutine ormas_mkden2x_old
!######################################################################
!######################################################################
subroutine ormas_mkden2x_old_bbp(nproc, cic, den2xp)

  use, intrinsic :: iso_c_binding
  use mod_const, only : czero, runit
  use mod_ormas, only : nact, det_allowed
  use mod_ormas, only : ndist_alph, nstr_alph, lstr_alph_dist, dist_str_alph, &
       & n1x_alph, p1x_alph, h1x_alph, eq1x_alph, sgn1x_alph
  use mod_ormas, only : ndist_beta, nstr_beta, lstr_beta_dist, dist_str_beta, &
       & n1x_beta, p1x_beta, h1x_beta, eq1x_beta, sgn1x_beta
  use mod_ormas, only : n1x_m_alph, map1x_m_alph
  use mod_ormas, only : n1x_m_beta, map1x_m_beta

  implicit none
  integer(c_long), intent(in) :: nproc
  complex(c_double_complex), intent(in) :: cic(1:nstr_alph, 1:nstr_beta)
  complex(c_double_complex), intent(inout) :: den2xp(1:nact, 1:nact, 1:nact, 1:nact, 0:(nproc-1))

  integer(c_long), external :: util_omp_iproc
  integer(c_long) :: nsafe, idist, jdist, kdist, ldist, istr, jstr, kstr, lstr, &
       & i1x, j1x, iact, jact, kact, lact, lla, ula, sgn, iproc, i1x_m, j1x_m
  complex(c_double_complex) :: tmp
  complex(c_double_complex), allocatable :: scic(:,:,:)

  nsafe = nact * nact
  allocate(scic(1:nstr_alph, 1:nsafe, 0:(nproc-1)))

  !$omp parallel default(shared) private(iact,jact,kact,lact,kstr,lstr,&
  !$omp & idist,kdist,ldist,sgn,tmp,lla,ula,iproc,i1x,i1x_m,j1x,j1x_m)
  iproc = util_omp_iproc()
  !$omp do
  do istr = 1, nstr_beta
     idist = dist_str_beta(1, istr)   
     scic(1:nstr_alph, 1:nsafe, iproc) = czero
!M-adapt
     do i1x = 1, n1x_beta(0, istr)
!     do i1x_m = 1, n1x_m_beta(0, istr)
!        i1x = map1x_m_beta(i1x_m, 0, istr)
!M-adapt
        sgn = sgn1x_beta(i1x, istr)
        kstr = eq1x_beta(i1x, istr)
        kdist = dist_str_beta(1, kstr)
        do jdist = 1, ndist_alph
           if (det_allowed(jdist, idist) == 0 .or. &
             & det_allowed(jdist, kdist) == 0) cycle

           lla = lstr_alph_dist(1, jdist)
           ula = lstr_alph_dist(2, jdist)
           do jstr = lla, ula
              scic(jstr, i1x, iproc) = cic(jstr, kstr) * sgn
           end do
        end do
     end do

!M-adapt
     do i1x = 1, n1x_beta(0, istr)
!     do i1x_m = 1, n1x_m_beta(0, istr)
!        i1x = map1x_m_beta(i1x_m, 0, istr)
!M-adapt
        iact = h1x_beta(i1x, istr) ! bra-hole
        jact = p1x_beta(i1x, istr) ! bra-particle
        kstr = eq1x_beta(i1x, istr)
        kdist = dist_str_beta(1, kstr)

!M-adapt
        do j1x = 1, n1x_beta(0, istr)
!        do j1x_m = 1, n1x_m_beta(0, istr)
!           j1x = map1x_m_beta(j1x_m, 0, istr)
!M-adapt
           kact = h1x_beta(j1x, istr) ! bra-hole
           lact = p1x_beta(j1x, istr) ! bra-particle
           lstr = eq1x_beta(j1x, istr)
           ldist = dist_str_beta(1, lstr)

           tmp = czero
           do jdist = 1, ndist_alph
              if (det_allowed(jdist, idist) == 0 .or. &
                & det_allowed(jdist, kdist) == 0 .or. &
                & det_allowed(jdist, ldist) == 0) cycle

              lla = lstr_alph_dist(1, jdist)
              ula = lstr_alph_dist(2, jdist)
              do jstr = lla, ula
                 tmp = tmp + conjg(scic(jstr, i1x, iproc)) &
                               & * scic(jstr, j1x, iproc)
              end do
           end do
           den2xp(iact, jact, kact, lact, iproc) = &   ! NOTE the sign
         & den2xp(iact, jact, kact, lact, iproc) + tmp ! NOTE the sign
        end do
     end do
  end do
!$omp end do
!$omp end parallel

  deallocate(scic)

end subroutine ormas_mkden2x_old_bbp
!######################################################################
subroutine ormas_mkden2x_old_abp(nproc, cic, den2xp)

  use, intrinsic :: iso_c_binding
  use mod_sph, only : mmax2
  use mod_const, only : czero, runit
  use mod_ormas, only : nact, det_allowed
  use mod_ormas, only : ndist_alph, nstr_alph, lstr_alph_dist, dist_str_alph, &
       & n1x_alph, p1x_alph, h1x_alph, eq1x_alph, sgn1x_alph
  use mod_ormas, only : ndist_beta, nstr_beta, lstr_beta_dist, dist_str_beta, &
       & n1x_beta, p1x_beta, h1x_beta, eq1x_beta, sgn1x_beta
  use mod_ormas, only : n1x_m_alph, map1x_m_alph
  use mod_ormas, only : n1x_m_beta, map1x_m_beta

  implicit none
  integer(c_long), intent(in) :: nproc
  complex(c_double_complex), intent(in) :: cic(1:nstr_alph, 1:nstr_beta)
  complex(c_double_complex), intent(inout) :: den2xp(1:nact, 1:nact, 1:nact, 1:nact, 0:(nproc-1))

  integer(c_long), external :: util_omp_iproc
  integer(c_long) :: istr, jstr, kstr, lstr, i1x, j1x, iact, jact, kact, lact, i1x_m, j1x_m
  integer(c_long) :: itype, lla, ula, llb, ulb, sgn1, sgn2, idist, jdist, kdist, iproc, m_ij, m_kl
  complex(c_double_complex) :: tmp

  do idist = 1, ndist_alph
     lla = lstr_alph_dist(1, idist)
     ula = lstr_alph_dist(2, idist)

!$omp parallel default(shared) private(iact,jact,kact,lact,kstr,lstr,llb,ulb,sgn1,sgn2,tmp,iproc,i1x,j1x,kdist,m_ij,m_kl)
     iproc = util_omp_iproc()
!$omp do
     do istr = lla, ula
!M-adapt
        do i1x = 1, n1x_alph(0, istr)
!        do m_ij = -mmax2, mmax2
!           m_kl = -m_ij
!        do i1x_m = 1, n1x_m_alph(m_ij, istr)
!           i1x = map1x_m_alph(i1x_m, m_ij, istr)
!M-adapt
           iact = h1x_alph (i1x, istr)
           jact = p1x_alph (i1x, istr)
           kstr = eq1x_alph(i1x, istr)
           sgn1 = sgn1x_alph(i1x, istr)
           kdist = dist_str_alph(1, kstr)
           do jdist = 1, ndist_beta
              if (det_allowed(idist, jdist) == 0) cycle
              if (det_allowed(kdist, jdist) == 0) cycle
              llb = lstr_beta_dist(1, jdist)
              ulb = lstr_beta_dist(2, jdist)
              do jstr = llb, ulb
                 tmp = conjg(cic(kstr, jstr)) * sgn1
!M-adapt
                 do j1x = 1, n1x_beta(0, jstr)
!                 do j1x_m = 1, n1x_m_beta(m_kl, jstr)
!                    j1x = map1x_m_beta(j1x_m, m_kl, jstr)
!M-adapt
                    kact = h1x_beta(j1x, jstr)
                    lact = p1x_beta(j1x, jstr)
                    lstr = eq1x_beta(j1x, jstr)
                    sgn2 = sgn1x_beta(j1x, jstr)
                    den2xp(iact, jact, kact, lact, iproc) = &                            ! NOTE the sign
                  & den2xp(iact, jact, kact, lact, iproc) + tmp * sgn2 * cic(istr, lstr) ! NOTE the sign
                 end do
              end do
           end do
        end do
!        end do
     end do
!$omp end do
!$omp end parallel
  end do

end subroutine ormas_mkden2x_old_abp
!######################################################################
