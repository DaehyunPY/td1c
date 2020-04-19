!######################################################################
subroutine ormas_spin_ras(cic, sz, s2)

  use, intrinsic :: iso_c_binding
  use mod_const, only : half
  use mod_ormas, only : nact, nelact

  implicit none
  complex(c_double_complex), intent(in) :: cic(1:*)
  real(c_double), intent(out) :: sz, s2

  sz = (nelact(1) - nelact(2)) * half
  s2 = sz * (sz + 1)
  if (nact > 0) call ormas_spin_ras_s2(cic, s2)

end subroutine ormas_spin_ras
!######################################################################
subroutine ormas_spin_ras_s2(cic, s2)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : mtot
  use mod_sph, only : mmax2
  use mod_const, only : zero, czero
  use mod_ormas, only : det_allowed, ndist_alph, nstr_alph, ndetx, mapr_detx, ntot_alph_beta, &
       & dist_str_alph, &
       & nstr_beta, dist_str_beta, n1x_beta, p1x_beta, h1x_beta, eq1x_beta, &
       & sgn1x_beta, n1x_m_alph, n1x_m_beta, map1x_m_alph, map1x_m_beta, h1x_alph, p1x_alph, eq1x_alph, sgn1x_alph
  use mod_ormas, only : mval_beta,nstr_dist_m_alph,llstr_dist_m_alph

  implicit none
  complex(c_double_complex), intent(in) :: cic(1:ndetx)
  real(c_double), intent(inout) :: s2

  real(c_double) :: s2p
  complex(c_double_complex) :: tmp
  integer(c_long) :: istr,jstr,kstr,lstr,ifun,jfun,kfun,lfun,i1x_m,i1x,j1x_m,j1x,m_ij,m_kl
  integer(c_long) :: idist,jdist,kdist,ldist,lla,ula,idet

  !$omp parallel default(shared) private(jstr,istr,m_kl,i1x,ifun,jfun,j1x,kfun,lfun,kstr,lstr,tmp,s2p) reduction(+:s2)
  s2p = zero
  !$omp do
  do idet = 1, ndetx
     jstr = mapr_detx(1,idet)
     istr = mapr_detx(2,idet)
     !##### DEBUG #####
     if (det_allowed(dist_str_alph(1,jstr),dist_str_beta(1,istr)) == 0) cycle
     !##### DEBUG #####
     do m_ij = -mmax2, mmax2
        do i1x_m = 1, n1x_m_beta(m_ij,istr)
           i1x = map1x_m_beta(i1x_m,m_ij,istr)
           ifun = h1x_beta (i1x,istr)
           jfun = p1x_beta (i1x,istr)
           lstr = eq1x_beta(i1x,istr)
           m_kl = -m_ij
           if (ifun == jfun) then
!              s2p = s2p + dble(conjg(cic(mapf_detx(jstr,istr))) &
!                                 & * cic(mapf_detx(jstr,lstr)))
!              s2p = s2p + dble(conjg(cic(mapf_detx(jstr,istr))) &
!                                 & * cic(mapf_detx(jstr,istr)))
              s2p = s2p + dble(conjg(cic(ntot_alph_beta(istr)+jstr)) &
                                 & * cic(ntot_alph_beta(istr)+jstr))
           end if

           tmp = czero
           do j1x_m = 1, n1x_m_alph(m_kl, jstr)
              j1x = map1x_m_alph(j1x_m, m_kl, jstr)
              kfun = h1x_alph (j1x, jstr)
              lfun = p1x_alph (j1x, jstr)
              kstr = eq1x_alph(j1x, jstr)
              if ((kfun /= jfun) .or. (lfun /= ifun)) cycle
              !### intensive if statement is not appropriate #####
              !### n1x_dist_m_aplh/beta should be developed ######
              if (det_allowed(dist_str_alph(1,kstr),dist_str_beta(1,lstr)) == 0) cycle
              !################################################
!              tmp = tmp + sgn1x_alph(j1x,jstr) * cic(mapf_detx(kstr,lstr))
              tmp = tmp + sgn1x_alph(j1x,jstr) * cic(ntot_alph_beta(lstr)+kstr)
           end do
!           s2p = s2p - dble(tmp) * sgn1x_beta(i1x,istr) * conjg(cic(mapf_detx(jstr,istr)))
!           s2p = s2p - dble(tmp * sgn1x_beta(i1x,istr) * conjg(cic(mapf_detx(jstr,istr))))
           s2p = s2p - dble(tmp * sgn1x_beta(i1x,istr) * conjg(cic(ntot_alph_beta(istr)+jstr)))
        end do
     end do
  end do
  !$omp end do
  s2 = s2 + s2p
  !$omp end parallel

end subroutine ormas_spin_ras_s2
!######################################################################
