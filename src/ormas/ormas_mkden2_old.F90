!######################################################################
subroutine ormas_mkden2_old(cic, den1, den2)

  use, intrinsic :: iso_c_binding
  use mod_const, only : czero, ctwo
  use mod_ormas, only : s2zero, nact, nstr_alph, nstr_beta

  implicit none
  complex(c_double_complex) , intent(in) :: cic(1:nstr_alph, 1:nstr_beta)
  complex(c_double_complex) , intent(in) :: den1(1:nact, 1:nact)
  complex(c_double_complex) , intent(out) :: den2(1:nact, 1:nact, 1:nact, 1:nact)

  integer(c_long) :: nact4, lwork
  integer(c_long) :: nproc,iproc,iact,jact,kact,lact
  complex(c_double_complex) , allocatable :: work(:,:)
  complex(c_double_complex) , allocatable :: den2p(:,:,:,:,:)
  integer(c_long), external :: util_omp_nproc

  nact4 = nact ** 4
  nproc = util_omp_nproc()
  lwork = max(nstr_alph, nstr_beta)
  allocate(work(lwork, 0:(nproc-1)))
  allocate(den2p(1:nact, 1:nact, 1:nact, 1:nact, 0:(nproc-1)))

  call util_zcopy(nact4, czero, 0, den2, 1)
  call util_zcopy(nact4*nproc, czero, 0, den2p, 1)

  call ormas_mkden2_old_bbp(cic, lwork, work, den2p, nproc)
  call ormas_mkden2_old_abp(cic, lwork, work, den2p, nproc)

  if (s2zero) then
     call util_zscal(nact4*nproc, ctwo, den2p, 1)
  else
     call ormas_mkden2_old_aap(cic, lwork, work, den2p, nproc)
     call ormas_mkden2_old_bap(cic, lwork, work, den2p, nproc)
  end if

  call ormas_mkden2_old_sum(den1, den2p, den2, nproc)

  deallocate(den2p)
  deallocate(work)

!debug
!  write(6, "('ormas_mkden2_old: den2')")
!  do lact = 1, nact
!  do kact = 1, nact
!  do jact = 1, nact
!  do iact = 1, nact
!     write(6, "(2f20.10)") den2(iact,jact,kact,lact)
!  end do
!  end do
!  end do
!  end do
!  stop
!debug

end subroutine ormas_mkden2_old
!######################################################################
!######################################################################
subroutine ormas_mkden2_old_bbp(cic, lwork, work, den2p, nproc)

  use, intrinsic :: iso_c_binding
  use mod_const, only : czero, runit, ctwo
  use mod_ormas, only : nact, det_allowed
  use mod_ormas, only : ndist_alph, nstr_alph, lstr_alph_dist
  use mod_ormas, only : nstr_beta, dist_str_beta, n1x_beta, p1x_beta, h1x_beta, eq1x_beta, sgn1x_beta

  implicit none
  integer(c_long), intent(in) :: lwork, nproc
  complex(c_double_complex) , intent(in) :: cic(1:nstr_alph, 1:nstr_beta)
  complex(c_double_complex) , intent(inout) :: work(1:lwork, 0:(nproc-1))
  complex(c_double_complex) , intent(inout) :: den2p(1:nact, 1:nact, 1:nact, 1:nact, 0:(nproc-1))

  integer(c_long) :: nact4, istr, jstr, kstr, lstr, i1x, j1x, iact, jact, kact, lact, &
       & iord, jord, idist, jdist, lla, ula, iproc
  integer(c_long), external :: util_omp_iproc

  nact4 = nact ** 4

!$omp parallel default(shared) private(iact,jact,kstr,iord,idist,lla,ula,kact,lact,lstr,jord,iproc)
     iproc = util_omp_iproc()
!$omp do
  do istr = 1, nstr_beta            ! loop over beta strings
     idist = dist_str_beta(1, istr)
     do i1x = 1, n1x_beta(0, istr)  ! loop over nonzero excitations within ormas space
        iact = h1x_beta (i1x, istr) ! bra-hole
        jact = p1x_beta (i1x, istr) ! bra-particle
        kstr = eq1x_beta(i1x, istr) ! matching string
        iord = nact * (iact - 1) + jact

        ! vectorize for alpha strings
        do jdist = 1, ndist_alph
!          if (det_allowed(jdist, idist) == 1) then
           if (det_allowed(jdist, idist) /= 0) then
              lla = lstr_alph_dist(1, jdist)
              ula = lstr_alph_dist(2, jdist)
              do jstr = lla, ula
                 work(jstr, iproc) = conjg(cic(jstr, istr)) * sgn1x_beta(i1x, istr)
              end do
           end if
        end do

        do j1x = 1, n1x_beta(0, kstr)  ! loop over excitations from kstr within ormas space
           kact = h1x_beta (j1x, kstr) ! bra-hole
           lact = p1x_beta (j1x, kstr) ! bra-particle
           lstr = eq1x_beta(j1x, kstr) ! matching string
           jord = nact * (kact - 1) + lact

           if (jord <= iord) then
              ! vectorize for alpha strings
              do jdist = 1, ndist_alph
!                if (det_allowed(jdist, idist) == 1) then
                 if (det_allowed(jdist, idist) /= 0) then
                    lla = lstr_alph_dist(1, jdist)
                    ula = lstr_alph_dist(2, jdist)
                    do jstr = lla, ula
                       den2p(jact, iact, lact, kact, iproc) = &
                     & den2p(jact, iact, lact, kact, iproc) + &
                     & work(jstr, iproc) * sgn1x_beta(j1x, kstr) * cic(jstr, lstr)
                    end do
                 end if
              end do
           end if
        end do
     end do
  end do
!$omp end do
!$omp end parallel

end subroutine ormas_mkden2_old_bbp
!######################################################################
!######################################################################
subroutine ormas_mkden2_old_aap(cic, lwork, work, den2p, nproc)

  use, intrinsic :: iso_c_binding
  use mod_const, only : czero, runit, ctwo
  use mod_ormas, only : nact, det_allowed
  use mod_ormas, only : ndist_alph, nstr_alph, lstr_alph_dist, &
       & n1x_alph, p1x_alph, h1x_alph, eq1x_alph, sgn1x_alph
  use mod_ormas, only : ndist_beta, nstr_beta, lstr_beta_dist

  implicit none
  integer(c_long), intent(in) :: lwork, nproc
  complex(c_double_complex) , intent(in) :: cic(1:nstr_alph, 1:nstr_beta)
  complex(c_double_complex) , intent(inout) :: work(1:lwork, 0:(nproc-1))
  complex(c_double_complex) , intent(inout) :: den2p(1:nact, 1:nact, 1:nact, 1:nact, 0:(nproc-1))

  integer(c_long) :: nact4, istr, jstr, kstr, lstr, i1x, j1x, iact, jact, kact, lact, &
       & iord, jord, idist, jdist, lla, ula, llb, ulb, iproc
  integer(c_long), external :: util_omp_iproc

  nact4 = nact ** 4

  do idist = 1, ndist_alph
     lla = lstr_alph_dist(1, idist)
     ula = lstr_alph_dist(2, idist)

!$omp parallel default(shared) private(istr,i1x,iact,jact,kstr,iord,jdist,llb,ulb,jstr,j1x,kact,lact,lstr,jord,iproc)
     iproc = util_omp_iproc()
!$omp do
     do istr = lla, ula                ! loop over alpha strings
        do i1x = 1, n1x_alph(0, istr)  ! loop over nonzero excitations from istra
           iact = h1x_alph (i1x, istr) ! bra-hole
           jact = p1x_alph (i1x, istr) ! bra-particle
           kstr = eq1x_alph(i1x, istr) ! matching string
           iord = nact * (iact - 1) + jact
  
           ! vectorize for beta strings
           do jdist = 1, ndist_beta
!             if (det_allowed(idist, jdist) == 1) then
              if (det_allowed(idist, jdist) /= 0) then
                 llb = lstr_beta_dist(1, jdist)
                 ulb = lstr_beta_dist(2, jdist)
                 do jstr = llb, ulb
                    work(jstr, iproc) = conjg(cic(istr, jstr)) * sgn1x_alph(i1x, istr)
                 end do
              end if
           end do
  
           do j1x = 1, n1x_alph(0, kstr)  ! loop over nonzero excitations from kstr
              kact = h1x_alph (j1x, kstr) ! bra-hole
              lact = p1x_alph (j1x, kstr) ! bra-particle
              lstr = eq1x_alph(j1x, kstr) ! matching string
              jord = nact * (kact - 1) + lact
  
              if (jord <= iord) then
                 ! vectorize for beta strings
                 do jdist = 1, ndist_beta
!                   if (det_allowed(idist, jdist) == 1) then
                    if (det_allowed(idist, jdist) /= 0) then
                       llb = lstr_beta_dist(1, jdist)
                       ulb = lstr_beta_dist(2, jdist)
                       do jstr = llb, ulb
                          den2p(jact, iact, lact, kact, iproc) = &
                        & den2p(jact, iact, lact, kact, iproc) + &
                        & work(jstr, iproc) * sgn1x_alph(j1x, kstr) * cic(lstr, jstr)
                       end do
                    end if
                 end do
              end if
           end do
        end do
     end do
!$omp end do
!$omp end parallel
  end do

end subroutine ormas_mkden2_old_aap
!######################################################################
!######################################################################
subroutine ormas_mkden2_old_abp(cic, lwork, work, den2p, nproc)

  use, intrinsic :: iso_c_binding
  use mod_const, only : czero, runit, ctwo
  use mod_ormas, only : nact, det_allowed
  use mod_ormas, only : ndist_alph, nstr_alph, lstr_alph_dist, &
       & n1x_alph, p1x_alph, h1x_alph, eq1x_alph, sgn1x_alph
  use mod_ormas, only : ndist_beta, nstr_beta, lstr_beta_dist, &
       & n1x_beta, p1x_beta, h1x_beta, eq1x_beta, sgn1x_beta

  implicit none
  integer(c_long), intent(in) :: lwork, nproc
  complex(c_double_complex) , intent(in) :: cic(1:nstr_alph, 1:nstr_beta)
  complex(c_double_complex) , intent(inout) :: work(1:lwork, 0:(nproc-1))
  complex(c_double_complex) , intent(inout) :: den2p(1:nact, 1:nact, 1:nact, 1:nact, 0:(nproc-1))

  integer(c_long) :: nact4, istr, jstr, kstr, lstr, i1x, j1x, iact, jact, kact, lact, &
       & iord, jord, idist, jdist, lla, ula, llb, ulb, iproc
  integer(c_long), external :: util_omp_iproc

  nact4 = nact ** 4

  do idist = 1, ndist_alph
     lla = lstr_alph_dist(1, idist)
     ula = lstr_alph_dist(2, idist)

!$omp parallel default(shared) private(iact,jact,kstr,iord,llb,ulb,kact,lact,lstr,jord,iproc)
     iproc = util_omp_iproc()
!$omp do
     do istr = lla, ula
        do i1x = 1, n1x_alph(0, istr)
           iact = h1x_alph (i1x, istr)
           jact = p1x_alph (i1x, istr)
           kstr = eq1x_alph(i1x, istr)
           iord = nact * (iact - 1) + jact

           do jdist = 1, ndist_beta
              if (det_allowed(idist, jdist) == 0) cycle
              llb = lstr_beta_dist(1, jdist)
              ulb = lstr_beta_dist(2, jdist)
              do jstr = llb, ulb
                 work(jstr, iproc) = conjg(cic(istr, jstr)) * sgn1x_alph(i1x, istr)
              end do

              do jstr = llb, ulb
                 do j1x = 1, n1x_beta(0, jstr)
                    kact = h1x_beta (j1x, jstr)
                    lact = p1x_beta (j1x, jstr)
                    lstr = eq1x_beta(j1x, jstr)
                    jord = nact * (kact - 1) + lact

                    if (jord <= iord) then
                       den2p(jact, iact, lact, kact, iproc) = &
                     & den2p(jact, iact, lact, kact, iproc) + &
                     & work(jstr, iproc) * sgn1x_beta(j1x, jstr) * cic(kstr, lstr)
                    end if
                 end do
              end do
           end do
        end do
     end do
!$omp end do
!$omp end parallel
  end do

end subroutine ormas_mkden2_old_abp
!######################################################################
!######################################################################
subroutine ormas_mkden2_old_bap(cic, lwork, work, den2p, nproc)

  use, intrinsic :: iso_c_binding
  use mod_const, only : czero, runit, ctwo
  use mod_ormas, only : nact, det_allowed
  use mod_ormas, only : ndist_alph, nstr_alph, lstr_alph_dist, &
       & n1x_alph, p1x_alph, h1x_alph, eq1x_alph, sgn1x_alph
  use mod_ormas, only : ndist_beta, nstr_beta, lstr_beta_dist, &
       & n1x_beta, p1x_beta, h1x_beta, eq1x_beta, sgn1x_beta

  implicit none
  integer(c_long), intent(in) :: lwork, nproc
  complex(c_double_complex) , intent(in) :: cic(1:nstr_alph, 1:nstr_beta)
  complex(c_double_complex) , intent(inout) :: work(1:lwork, 0:(nproc-1))
  complex(c_double_complex) , intent(inout) :: den2p(1:nact, 1:nact, 1:nact, 1:nact, 0:(nproc-1))

  integer(c_long) :: nact4, istr, jstr, kstr, lstr, i1x, j1x, iact, jact, kact, lact, &
       & iord, jord, idist, jdist, lla, ula, llb, ulb, iproc
  integer(c_long), external :: util_omp_iproc

  nact4 = nact ** 4

  do idist = 1, ndist_beta
     llb = lstr_beta_dist(1, idist)
     ulb = lstr_beta_dist(2, idist)

!$omp parallel default(shared) private(istr,i1x,iact,jact,kstr,iord,jdist,lla,ula,jstr,j1x,kact,lact,lstr,jord,iproc)
     iproc = util_omp_iproc()
!$omp do
     do istr = llb, ulb
        do i1x = 1, n1x_beta(0, istr)
           iact = h1x_beta (i1x, istr)
           jact = p1x_beta (i1x, istr)
           kstr = eq1x_beta(i1x, istr)
           iord = nact * (iact - 1) + jact     

           do jdist = 1, ndist_alph
!             if (det_allowed(jdist, idist) == 1) then
              if (det_allowed(jdist, idist) /= 0) then
                 lla = lstr_alph_dist(1, jdist)
                 ula = lstr_alph_dist(2, jdist)
                 do jstr = lla, ula
                    work(jstr, iproc) = conjg(cic(jstr, istr)) * sgn1x_beta(i1x, istr)
                 end do
  
                 do jstr = lla, ula
                    do j1x = 1, n1x_alph(0, jstr)
                       kact = h1x_alph (j1x, jstr)
                       lact = p1x_alph (j1x, jstr)
                       lstr = eq1x_alph(j1x, jstr)
                       jord = nact * (kact - 1) + lact
              
                       if (jord <= iord) then     
                          den2p(jact, iact, lact, kact, iproc) = &
                        & den2p(jact, iact, lact, kact, iproc) + &
                        & work(jstr, iproc) * sgn1x_alph(j1x, jstr) * cic(lstr, kstr)
                       end if
                    end do
                 end do
              end if
           end do
        end do
     end do
!$omp end do
!$omp end parallel
  end do

end subroutine ormas_mkden2_old_bap
!######################################################################
!######################################################################
subroutine ormas_mkden2_old_sum(den1, den2p, den2, nproc)

  use, intrinsic :: iso_c_binding
  use mod_const, only : runit
  use mod_ormas, only : nact

  implicit none
  integer(c_long), intent(in) :: nproc
  complex(c_double_complex) , intent(in) :: den1(1:nact, 1:nact)
  complex(c_double_complex) , intent(in) :: den2p(1:nact, 1:nact, 1:nact, 1:nact, 0:(nproc-1))
  complex(c_double_complex) , intent(out) :: den2(1:nact, 1:nact, 1:nact, 1:nact)

  integer(c_long) :: nact4, iact, jact, kact, lact, iord, jord, iproc

  nact4 = nact ** 4
  do iproc = 0, nproc - 1
     call zaxpy(nact4, runit, den2p(1,1,1,1,iproc), 1, den2, 1)
!     call util_zaxpy(nact4, runit, den2p(1,1,1,1,iproc), 1, den2, 1)
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
              if (jord > iord) then
                 den2(jact, iact, lact, kact) = &
               & den2(lact, kact, jact, iact)
              end if
           end do
        end do
     end do
  end do

end subroutine ormas_mkden2_old_sum
!######################################################################
