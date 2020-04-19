!######################################################################
subroutine ormas_spin_old(cic, sz, s2)

  use, intrinsic :: iso_c_binding
  use mod_const, only : half
  use mod_ormas, only : nelact

  implicit none
  complex(c_double_complex), intent(in) :: cic(1:*)
  real(c_double), intent(out) :: sz, s2

  sz = (nelact(1) - nelact(2)) * half
  s2 = sz * (sz + 1)
  call ormas_spin_old_s2(cic, s2)

end subroutine ormas_spin_old
!######################################################################
subroutine ormas_spin_old_s2(cic, s2)

  use, intrinsic :: iso_c_binding
  use mod_const, only : zero, runit
  use mod_ormas, only : det_allowed, ndist_alph, nstr_alph, lstr_alph_dist, &
       & dist_str_alph, nstr_beta, dist_str_beta, n1x_beta, p1x_beta, h1x_beta, &
       & eq1x_beta, sgn1x_beta
!OLD  use mod_ormas, only : n1xr_alph, l1xr_alph, r1xr_alph, sgn1xr_alph, 

  implicit none
  complex(c_double_complex), intent(in) :: cic(1:nstr_alph, 1:nstr_beta)
  real(c_double), intent(inout) :: s2

  stop "ormas_spin_old_s2 no longer supported."

!OLD  real(c_double) :: s2p
!OLD  complex(c_double_complex) :: sgnb, sgna
!OLD  integer(c_long) :: istr, jstr, kstr, lstr, i1x, ifun, jfun, &
!OLD       & idist, jdist, ii, n1xra, lla, ula, iproc
!OLD  integer(c_long), external :: util_omp_iproc
!OLD
!OLD!$omp parallel default(shared) private(istr,jstr,kstr,lstr,ii, &
!OLD!$omp & idist,jdist,i1x,ifun,jfun,sgnb,sgna,lla,ula,s2p,n1xra,iproc) &
!OLD!$omp & reduction(+:s2)
!OLD  iproc = util_omp_iproc()
!OLD  s2p = zero
!OLD!$omp do
!OLD  do istr = 1, nstr_beta
!OLD     idist = dist_str_beta(1, istr)
!OLD     do i1x = 1, n1x_beta(0, istr)
!OLD        ifun = h1x_beta  (i1x, istr)
!OLD        jfun = p1x_beta  (i1x, istr)
!OLD        kstr = eq1x_beta (i1x, istr)
!OLD        sgnb = sgn1x_beta(i1x, istr) * runit
!OLD
!OLD        if (ifun == jfun) then
!OLD           do jdist = 1, ndist_alph
!OLD              if (det_allowed(jdist, idist) == 0) cycle
!OLD              lla = lstr_alph_dist(1, jdist)
!OLD              ula = lstr_alph_dist(2, jdist)
!OLD              do jstr = lla, ula
!OLD                 s2p = s2p + dble(conjg(cic(jstr, istr)) &
!OLD                                    & * cic(jstr, istr))
!OLD              end do
!OLD           end do
!OLD        end if
!OLD
!OLD        n1xra = n1xr_alph(jfun, ifun)
!OLD        do ii = 1, n1xra
!OLD           jstr = r1xr_alph(ii, jfun, ifun)
!OLD           jdist = dist_str_alph(1, jstr)
!OLD           if (det_allowed(jdist, idist) == 0) cycle
!OLD
!OLD           lstr = l1xr_alph(ii, jfun, ifun)
!OLD           sgna = sgn1xr_alph(ii, jfun, ifun) * runit
!OLD           s2p = s2p - dble(conjg(cic(jstr, istr)) &
!OLD                              & * cic(lstr, kstr) * sgnb * sgna)
!OLD        end do
!OLD     end do
!OLD  end do
!OLD!$omp end do
!OLD  s2 = s2 + s2p
!OLD!$omp end parallel

end subroutine ormas_spin_old_s2
!######################################################################
