!######################################################################
subroutine ormas_mkden1_old(cic, den1)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : smul
  use mod_const, only : czero, runit, ctwo
  use mod_ormas, only : nelact, nact, nstr_alph, nstr_beta

  implicit none
  complex(c_double_complex) , intent(in) :: cic(1:nstr_alph, 1:nstr_beta)
  complex(c_double_complex) , intent(out) :: den1(1:nact, 1:nact)

  integer(c_long) :: nact2, iproc, nproc, iact, jact
  complex(c_double_complex) , allocatable :: den1p(:,:,:)
  integer(c_long), external :: util_omp_nproc

  nact2 = nact * nact
  nproc = util_omp_nproc()
  allocate(den1p(1:nact, 1:nact, 0:(nproc-1)))
  call util_zcopy(nact2*nproc, czero, 0, den1p, 1)

  call ormas_mkden1_old_bb(cic, den1p, nproc)
  if (smul == 1 .and. nelact(1) == nelact(2)) then
!     call zscal(nact2*nproc, ctwo, den1p, 1)
     call util_zscal(nact2*nproc, ctwo, den1p, 1)
  else
     stop 'ormas_mkden1_old: only for NA = NB.'
!     call ormas_mkden1_old_aa(cic, den1p)
  end if

  call util_zcopy(nact2, czero, 0, den1, 1)
  do iproc = 0, nproc - 1
     call zaxpy(nact2, runit, den1p(1,1,iproc), 1, den1, 1)
!     call util_zaxpy(nact2, runit, den1p(1,1,iproc), 1, den1, 1)
  end do

  deallocate(den1p)

!debug
!  write(6, "('ormas_mkden1_old:')")
!  do iact = 1, nact
!     do jact = 1, nact
!        write(6, "(f20.10)", advance = 'no') dble(den1(jact, iact))
!     end do
!     write(6, *)
!  end do
!  stop
!debug
end subroutine ormas_mkden1_old
!######################################################################
!######################################################################
subroutine ormas_mkden1_old_bb(cic, den1p, nproc)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nact, det_allowed, ndist_alph
  use mod_ormas, only : lstr_alph_dist, dist_str_beta, nstr_alph
  use mod_ormas, only : nstr_beta, n1x_beta, p1x_beta, h1x_beta, eq1x_beta, sgn1x_beta
  use mod_ormas, only : n1x_m_alph, map1x_m_alph
  use mod_ormas, only : n1x_m_beta, map1x_m_beta

  implicit none
  integer(c_long), intent(in) :: nproc
  complex(c_double_complex) , intent(in) :: cic(1:nstr_alph, 1:nstr_beta)
  complex(c_double_complex) , intent(inout) :: den1p(1:nact, 1:nact, 0:(nproc-1))

  integer(c_long) :: istr, jstr, kstr, i1x, i1x_m, iact, jact, iproc
  integer(c_long) :: idist, jdist, kdist, lla, ula
  integer(c_long), external :: util_omp_iproc

!$omp parallel default(shared) private(istr,idist,i1x,i1x_m,iact,jact,kstr,kdist,jstr,jdist,lla,ula,iproc)
  iproc = util_omp_iproc()
!$omp do 
  do istr = 1, nstr_beta
     idist = dist_str_beta(1, istr)
!M-adapt
!    do i1x = 1, n1x_beta(0, istr)
     do i1x_m = 1, n1x_m_beta(0, istr)
        i1x = map1x_m_beta(i1x_m, 0, istr)
!M-adapt
        iact = h1x_beta (i1x, istr)
        jact = p1x_beta (i1x, istr)
        kstr = eq1x_beta(i1x, istr)
        kdist = dist_str_beta(1, kstr)
        do jdist = 1, ndist_alph
           if (det_allowed(jdist, idist) /= 0 .and. &
             & det_allowed(jdist, kdist) /= 0) then ! DIFFER FROM AAP
              lla = lstr_alph_dist(1, jdist)
              ula = lstr_alph_dist(2, jdist)
              do jstr = lla, ula
                 den1p(jact, iact, iproc) = den1p(jact, iact, iproc) &
                      & + conjg(cic(jstr, istr)) &
                      &       * cic(jstr, kstr) * sgn1x_beta(i1x, istr)
              end do
           end if
        end do
     end do
  end do
!$omp end do
!$omp end parallel

end subroutine ormas_mkden1_old_bb
!######################################################################
