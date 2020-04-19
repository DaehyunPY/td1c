!################################################################################
subroutine ormas_madapt(mtot)

  use, intrinsic :: iso_c_binding
  use mod_sph, only : mmax2
  use mod_ormas, only : iprint, nact, nfun, cic_old, mval
  use mod_ormas, only : nstr_alph, n1x_alph, p1x_alph, h1x_alph, &
       n1x_m_alph, map1x_m_alph, max1x_alph, max1x_m_alph
  use mod_ormas, only : nstr_beta, n1x_beta, p1x_beta, h1x_beta, &
       n1x_m_beta, map1x_m_beta, max1x_beta, max1x_m_beta
  
  implicit none
  integer(c_int), intent(in) :: mtot

!tmp  integer(c_int) :: max1x_alph, max1x_m_alph
!tmp  integer(c_int) :: max1x_beta, max1x_m_beta

  call ormas_max1x_m_spin(nstr_alph, max1x_alph, n1x_alph, p1x_alph, h1x_alph, max1x_m_alph)
  allocate(n1x_m_alph(-mmax2:mmax2, 1:nstr_alph))
  allocate(map1x_m_alph(1:max1x_m_alph, -mmax2:mmax2, 1:nstr_alph))
  call ormas_madapt_int1x_spin(nstr_alph, max1x_alph, n1x_alph, p1x_alph, h1x_alph, &
       max1x_m_alph, n1x_m_alph, map1x_m_alph)

  call ormas_max1x_m_spin(nstr_beta, max1x_beta, n1x_beta, p1x_beta, h1x_beta, max1x_m_beta)
  allocate(n1x_m_beta(-mmax2:mmax2, 1:nstr_beta))
  allocate(map1x_m_beta(1:max1x_m_beta, -mmax2:mmax2, 1:nstr_beta))
  call ormas_madapt_int1x_spin(nstr_beta, max1x_beta, n1x_beta, p1x_beta, h1x_beta, &
       max1x_m_beta, n1x_m_beta, map1x_m_beta)

  call ormas_madapt_str(mtot)

end subroutine ormas_madapt
!################################################################################
!################################################################################
subroutine ormas_max1x_m_spin(nstr, max1x, n1x, p1x, h1x, max1x_m)

  use, intrinsic :: iso_c_binding
  use mod_sph, only : mmax2
  use mod_ormas, only : ncore, nact, nfun, mval

  implicit none
  integer(c_int), intent(in) :: nstr
  integer(c_int), intent(in) :: max1x, n1x(-3:1, 1:*)
  integer(c_int), intent(in) :: p1x(1:nact*nact, 1:*)
  integer(c_int), intent(in) :: h1x(1:nact*nact, 1:*)
  integer(c_int), intent(out) :: max1x_m

  integer(c_int) :: istr, i1x, ifun, jfun, m_ij, i1x_m, n1xm(-mmax2:mmax2)

  max1x_m = 0
  do istr = 1, nstr
     n1xm(-mmax2:mmax2) = 0
     do i1x = 1, n1x(0,istr)
        ifun = p1x(i1x,istr) ! ket particle
        jfun = h1x(i1x,istr) ! ket hole
        m_ij = - mval(ncore+ifun) + mval(ncore+jfun)
        n1xm(m_ij) = n1xm(m_ij) + 1
     end do
     do m_ij = -mmax2, mmax2
        max1x_m = max(max1x_m, n1xm(m_ij))
     end do
  end do

!debug
!  write(6, "('max1x_m: ', 2i10)") max1x_m, nact*nact
!  stop
!debug

end subroutine ormas_max1x_m_spin
!################################################################################
subroutine ormas_madapt_int1x_spin(nstr, max1x, n1x, p1x, h1x, max1x_m, n1x_m, map1x_m)

  use, intrinsic :: iso_c_binding
  use mod_sph, only : mmax2
  use mod_ormas, only : ncore, nact, nfun, mval

  implicit none
  integer(c_int), intent(in) :: nstr
  integer(c_int), intent(in) :: max1x, max1x_m
  integer(c_int), intent(in) :: n1x(-3:1, 1:*)
  integer(c_int), intent(in) :: p1x(1:nact*nact, 1:*)
  integer(c_int), intent(in) :: h1x(1:nact*nact, 1:*)
  integer(c_int), intent(out) :: n1x_m(-mmax2:mmax2, 1:*)
  integer(c_int), intent(out) :: map1x_m(1:max1x_m, -mmax2:mmax2, 1:*)

  integer(c_int) :: istr, i1x, ifun, jfun, m_ij, i1x_m, maxn1x_m

  n1x_m(-mmax2:mmax2, 1:nstr) = 0
  map1x_m(1:max1x_m, -mmax2:mmax2, 1:nstr) = 0

  do istr = 1, nstr
     do i1x = 1, n1x(0,istr)
        ifun = p1x(i1x,istr) ! ket particle
        jfun = h1x(i1x,istr) ! ket hole
!BUG    m_ij = + mval(ncore + ifun) - mval(ncore + jfun)
        m_ij = - mval(ncore + ifun) + mval(ncore + jfun)
        n1x_m(m_ij, istr) = n1x_m(m_ij, istr) + 1
        map1x_m(n1x_m(m_ij, istr), m_ij, istr) = i1x
     end do
  end do

!debug     do istr = 1, nstr
!debug        do m_ij = -mmax2, mmax2
!debug           write(6, "(2i5)") istr, m_ij
!debug           do i1x_m = 1, n1x_m(m_ij, istr)
!debug              i1x = map1x_m(i1x_m, m_ij, istr)
!debug              ifun = p1x(i1x, istr)
!debug              jfun = h1x(i1x, istr)
!debug              write(6, "(5x, 4i5)") i1x_m, i1x, ifun, jfun
!debug           end do
!debug        end do
!debug     end do

end subroutine ormas_madapt_int1x_spin
!################################################################################
!################################################################################
subroutine ormas_madapt_str(mtot)

  use, intrinsic :: iso_c_binding
  use mod_sph, only : mmax2
  use mod_ormas, only : cic_old
  use mod_ormas, only : nstr_alph, orb_alph, nstr_alph_beta, llstr_alph_beta, ntot_alph_beta
  use mod_ormas, only : nstr_beta, orb_beta, nstr_beta_alph, llstr_beta_alph, ntot_beta_alph
  use mod_ormas, only : iprint, ncore, nact, nfun, mval, nelact, nsub, det_allowed
  use mod_ormas, only : ndist_alph, dist_str_alph, substr_alph, onv_alph, orb_alph
  use mod_ormas, only : ndist_beta, dist_str_beta, substr_beta, onv_beta, orb_beta
  use mod_ormas, only : mmin_alph, mmax_alph, mval_alph
  use mod_ormas, only : mmin_beta, mmax_beta, mval_beta
  use mod_ormas, only : map2to1_alph, map1to2_alph
  use mod_ormas, only : map2to1_beta, map1to2_beta
  use mod_ormas, only : ndet, ndetx, mapr_detx
  use mod_ormas, only : n1x_alph, p1x_alph, h1x_alph, eq1x_alph, sgn1x_alph, max1x_m_alph, n1x_m_alph, map1x_m_alph
  use mod_ormas, only : n1x_beta, p1x_beta, h1x_beta, eq1x_beta, sgn1x_beta, max1x_m_beta, n1x_m_beta, map1x_m_beta
!OLD  use mod_ormas, only : n1xr_alph, r1xr_alph, l1xr_alph, sgn1xr_alph
!OLD  use mod_ormas, only : n1xr_beta, r1xr_beta, l1xr_beta, sgn1xr_beta
  use mod_ormas, only : nstr_alph_dist, lstr_alph_dist
  use mod_ormas, only : nstr_beta_dist, lstr_beta_dist
  use mod_ormas, only : nstr_dist_m_alph, llstr_dist_m_alph
  use mod_ormas, only : nstr_dist_m_beta, llstr_dist_m_beta

  implicit none
  integer(c_int), intent(in) :: mtot
  logical(c_bool) :: found
  integer(c_int) :: isub, istr, jstr, ifun, jfun, iel, idet, ndety, mvalc, tmp, ii, iact, jact
  integer(c_int) :: mval_str, istr_m, istr2, jstr2, idist, jdist, i1x, i1x_m, m_ij, lls, uls, mstra, mstrb
  integer(c_int) :: max_nstr_m_alph, max_nstr_m_beta
  integer(c_int), allocatable :: tocc(:), mval_alph2(:), mval_beta2(:)
  integer(c_int), allocatable :: strmap_m_alph(:,:), str_dist_m_alph(:,:,:)
  integer(c_int), allocatable :: strmap_m_beta(:,:), str_dist_m_beta(:,:,:)
  integer(c_int), allocatable :: dist_str_alph2(:,:), onv_alph2(:,:), orb_alph2(:,:)
  integer(c_int), allocatable :: dist_str_beta2(:,:), onv_beta2(:,:), orb_beta2(:,:)
  integer(c_int), allocatable :: n1x_alph2(:,:), p1x_alph2(:,:), h1x_alph2(:,:)
  integer(c_int), allocatable :: eq1x_alph2(:,:), sgn1x_alph2(:,:), n1x_m_alph2(:,:), map1x_m_alph2(:,:,:)
  integer(c_int), allocatable :: n1x_beta2(:,:), p1x_beta2(:,:), h1x_beta2(:,:)
  integer(c_int), allocatable :: eq1x_beta2(:,:), sgn1x_beta2(:,:), n1x_m_beta2(:,:), map1x_m_beta2(:,:,:)
!OLD  integer(c_int), allocatable :: n1xr_alph2(:,:), r1xr_alph2(:,:,:), l1xr_alph2(:,:,:), sgn1xr_alph2(:,:,:)
!OLD  integer(c_int), allocatable :: n1xr_beta2(:,:), r1xr_beta2(:,:,:), l1xr_beta2(:,:,:), sgn1xr_beta2(:,:,:)
  

!debug
  integer(c_int) :: it1,it2,it3,lla,ula
!debug


  allocate(nstr_alph_beta(1:nstr_beta))
  allocate(nstr_beta_alph(1:nstr_alph))
  allocate(llstr_alph_beta(1:nstr_beta))
  allocate(llstr_beta_alph(1:nstr_alph))
  allocate(ntot_alph_beta(1:nstr_beta))
  allocate(ntot_beta_alph(1:nstr_alph))
  allocate(mval_alph(1:nstr_alph))
  allocate(mval_beta(1:nstr_beta))
  allocate(map2to1_alph(1:nstr_alph))
  allocate(map2to1_beta(1:nstr_beta))
  allocate(map1to2_alph(1:nstr_alph))
  allocate(map1to2_beta(1:nstr_beta))

  ! scratch arrays for reordering
  allocate(mval_alph2(1:nstr_alph))
  allocate(mval_beta2(1:nstr_beta))
  allocate(dist_str_alph2(1:2, 1:nstr_alph))
  allocate(dist_str_beta2(1:2, 1:nstr_beta))
  allocate(onv_alph2(1:nact, 1:nstr_alph))
  allocate(onv_beta2(1:nact, 1:nstr_beta))
  allocate(orb_alph2(0:nelact(1), 1:nstr_alph))
  allocate(orb_beta2(0:nelact(2), 1:nstr_beta))
  allocate(n1x_alph2(-3:1, 1:nstr_alph))
  allocate(p1x_alph2(1:nact*nact, 1:nstr_alph))
  allocate(h1x_alph2(1:nact*nact, 1:nstr_alph))
  allocate(eq1x_alph2(1:nact*nact, 1:nstr_alph))
  allocate(sgn1x_alph2(1:nact*nact, 1:nstr_alph))
  allocate(n1x_m_alph2(-mmax2:mmax2, 1:nstr_alph))
  allocate(map1x_m_alph2(1:max1x_m_alph, -mmax2:mmax2, 1:nstr_alph))
  allocate(n1x_beta2(-3:1, 1:nstr_beta))
  allocate(p1x_beta2(1:nact*nact, 1:nstr_beta))
  allocate(h1x_beta2(1:nact*nact, 1:nstr_beta))
  allocate(eq1x_beta2(1:nact*nact, 1:nstr_beta))
  allocate(sgn1x_beta2(1:nact*nact, 1:nstr_beta))
  allocate(n1x_m_beta2(-mmax2:mmax2, 1:nstr_beta))
  allocate(map1x_m_beta2(1:max1x_m_beta, -mmax2:mmax2, 1:nstr_beta))
!OLD  allocate(n1xr_alph2(             1:nact, 1:nact))
!OLD  allocate(r1xr_alph2(1:nstr_alph, 1:nact, 1:nact))
!OLD  allocate(l1xr_alph2(1:nstr_alph, 1:nact, 1:nact))
!OLD  allocate(sgn1xr_alph2(1:nstr_alph, 1:nact, 1:nact))
!OLD  allocate(n1xr_beta2(             1:nact, 1:nact))
!OLD  allocate(r1xr_beta2(1:nstr_beta, 1:nact, 1:nact))
!OLD  allocate(l1xr_beta2(1:nstr_beta, 1:nact, 1:nact))
!OLD  allocate(sgn1xr_beta2(1:nstr_beta, 1:nact, 1:nact))

  !##### total M of each spin string
  mvalc = 0
  do ifun = 1, ncore
     mvalc = mvalc + mval(ifun)
  end do
  mmin_alph = 0
  mmax_alph = 0
  do istr = 1, nstr_alph
     tmp = 0
     do iel = 1, nelact(1)
        tmp = tmp + mval(ncore+orb_alph(iel,istr))
     end do
     mval_alph(istr) = mvalc + tmp
     if (mval_alph(istr) < mmin_alph) mmin_alph = mval_alph(istr)
     if (mval_alph(istr) > mmax_alph) mmax_alph = mval_alph(istr)
  end do
  mmin_beta = 0
  mmax_beta = 0
  do istr = 1, nstr_beta
     tmp = 0
     do iel = 1, nelact(2)
        tmp = tmp + mval(ncore+orb_beta(iel,istr))
     end do
     mval_beta(istr) = mvalc + tmp
     if (mval_beta(istr) < mmin_beta) mmin_beta = mval_beta(istr)
     if (mval_beta(istr) > mmax_beta) mmax_beta = mval_beta(istr)
  end do
  !DEBUG
  !write(6, "('boundaries of M_det:')")
  !write(6, "(2i5)") mmin_alph,mmin_beta
  !write(6, "(2i5)") mmax_alph,mmax_beta
  !stop
  !DEBUG

  ! number of strings for each distribution for each M
  allocate(nstr_dist_m_alph(1:ndist_alph, mmin_alph:mmax_alph))
  allocate(nstr_dist_m_beta(1:ndist_beta, mmin_beta:mmax_beta))
  allocate(str_dist_m_alph(1:nstr_alph, 1:ndist_alph, mmin_alph:mmax_alph)) ! temporary
  allocate(str_dist_m_beta(1:nstr_beta, 1:ndist_beta, mmin_beta:mmax_beta)) ! temporary

  nstr_dist_m_alph(1:ndist_alph, mmin_alph:mmax_alph) = 0
  do istr = 1, nstr_alph
     nstr_dist_m_alph(dist_str_alph(1,istr), mval_alph(istr)) = &
     nstr_dist_m_alph(dist_str_alph(1,istr), mval_alph(istr)) + 1
     str_dist_m_alph(nstr_dist_m_alph(dist_str_alph(1,istr),mval_alph(istr)),dist_str_alph(1,istr),mval_alph(istr)) = istr
  end do
  nstr_dist_m_beta(1:ndist_beta, mmin_beta:mmax_beta) = 0
  do istr = 1, nstr_beta
     nstr_dist_m_beta(dist_str_beta(1, istr), mval_beta(istr)) = &
     nstr_dist_m_beta(dist_str_beta(1, istr), mval_beta(istr)) + 1
     str_dist_m_beta(nstr_dist_m_beta(dist_str_beta(1,istr),mval_beta(istr)),dist_str_beta(1,istr),mval_beta(istr)) = istr
  end do
  !DEBUG
  !write(6, "('M- and distribution-resolved alpha strings:')")
  !do mval_str = mmin_alph, mmax_alph
  !   do idist = 1, ndist_alph
  !      do ii = 1, nstr_dist_m_alph(idist,mval_str)
  !         write(6, "(4i5)") mval_str,idist,ii,str_dist_m_alph(ii,idist,mval_str)
  !      end do
  !   end do
  !end do
  !write(6, "('M- and distribution-resolved beta strings:')")
  !do mval_str = mmin_beta, mmax_beta
  !   do idist = 1, ndist_beta
  !      do ii = 1, nstr_dist_m_beta(idist,mval_str)
  !         write(6, "(4i5)") mval_str,idist,ii,str_dist_m_beta(ii,idist,mval_str)
  !      end do
  !   end do
  !end do
  !stop
  !DEBUG

  ! map for M-ordered string to distribution-ordered string
  istr2 = 0
  do mval_str = mmin_alph, mmax_alph
     do idist = 1, ndist_alph
        do ii = 1, nstr_dist_m_alph(idist,mval_str)
           istr2 = istr2 + 1
           istr = str_dist_m_alph(ii,idist,mval_str)
           map2to1_alph(istr2) = istr
           map1to2_alph(istr) = istr2
        end do
     end do
  end do
  if (istr2 .ne. nstr_alph) stop "algorithum failure in ormas_madapt_str (1)."
  istr2 = 0
  do mval_str = mmin_beta, mmax_beta
     do idist = 1, ndist_beta
        do ii = 1, nstr_dist_m_beta(idist,mval_str)
           istr2 = istr2 + 1
           istr = str_dist_m_beta(ii,idist,mval_str)
           map2to1_beta(istr2) = istr
           map1to2_beta(istr) = istr2
        end do
     end do
  end do
  if (istr2 .ne. nstr_beta) stop "algorithum failure in ormas_madapt_str (2)."
  !DEBUG
  !write(6, "('alpha reordered strings:')")
  !do istr2 = 1, nstr_alph
  !   istr = map2to1_alph(istr2)
  !   idist = dist_str_alph(1,istr)
  !   mval_str = mval_alph(istr)
  !   write(6, "(4i10,2x)", advance = 'no') istr2,istr,idist,mval_str
  !   call ormas_occvec_print(6,.true.,nact,onv_alph(1,istr))
  !end do
  !write(6, "('beta reordered strings:')")
  !do istr2 = 1, nstr_beta
  !   istr = map2to1_beta(istr2)
  !   idist = dist_str_beta(1,istr)
  !   mval_str = mval_beta(istr)
  !   write(6, "(4i10,2x)", advance = 'no') istr2,istr,idist,mval_str
  !   call ormas_occvec_print(6,.true.,nact,onv_beta(1,istr))
  !end do
  !stop
  !DEBUG
  deallocate(str_dist_m_alph) ! temporary
  deallocate(str_dist_m_beta) ! temporary

  ! left boundary of strings for each distribution for each M
  allocate(llstr_dist_m_alph(1:ndist_alph, mmin_alph:mmax_alph))
  allocate(llstr_dist_m_beta(1:ndist_beta, mmin_beta:mmax_beta))
  tmp = 0
  do mval_str = mmin_alph, mmax_alph
     do idist = 1, ndist_alph
        llstr_dist_m_alph(idist,mval_str) = tmp + 1
        tmp = tmp + nstr_dist_m_alph(idist,mval_str)
     end do
  end do
  tmp = 0
  do mval_str = mmin_beta, mmax_beta
     do idist = 1, ndist_beta
        llstr_dist_m_beta(idist,mval_str) = tmp + 1
        tmp = tmp + nstr_dist_m_beta(idist,mval_str)
     end do
  end do

!  write(6, "('# ORMAS: test 1')")
  ! Now various pointers are reordered

  p1x_alph2 = 0
  h1x_alph2 = 0
  eq1x_alph2 = 0
  sgn1x_alph2 = 0
  do istr2 = 1, nstr_alph
     istr = map2to1_alph(istr2)
     mval_alph2(istr2) = mval_alph(istr)
     dist_str_alph2(1,istr2) = dist_str_alph(1,istr)

     onv_alph2(1:nact, istr2) = onv_alph(1:nact, istr)
     orb_alph2(0:nelact(1), istr2) = orb_alph(0:nelact(1), istr)

     n1x_alph2(-3:1,istr2) = n1x_alph(-3:1,istr)
     do i1x = 1, n1x_alph2(1,istr2)
        p1x_alph2(i1x,istr2) = p1x_alph(i1x,istr)
        h1x_alph2(i1x,istr2) = h1x_alph(i1x,istr)
        sgn1x_alph2(i1x,istr2) = sgn1x_alph(i1x,istr)
        if (eq1x_alph(i1x,istr) > 0) then
           eq1x_alph2(i1x,istr2) = map1to2_alph(eq1x_alph(i1x,istr))
        else
           eq1x_alph2(i1x,istr2) = eq1x_alph(i1x,istr)
        end if
     end do
     n1x_m_alph2(-mmax2:mmax2, istr2) = n1x_m_alph(-mmax2:mmax2, istr)
     map1x_m_alph2(1:max1x_m_alph, -mmax2:mmax2, istr2) = map1x_m_alph(1:max1x_m_alph, -mmax2:mmax2, istr)
  end do

!OLD  n1xr_alph2(1:nact,1:nact) = n1xr_alph(1:nact,1:nact)
!OLD  r1xr_alph2 = 0
!OLD  l1xr_alph2 = 0
!OLD  sgn1xr_alph2 = 0
!OLD  do iact = 1, nact
!OLD  do jact = 1, nact
!OLD     do ii = 1, n1xr_alph2(jact,iact)
!OLD        r1xr_alph2(ii,jact,iact) = map1to2_alph(r1xr_alph(ii,jact,iact))
!OLD        l1xr_alph2(ii,jact,iact) = map1to2_alph(l1xr_alph(ii,jact,iact))
!OLD        sgn1xr_alph2(ii,jact,iact) = sgn1xr_alph(ii,jact,iact)
!OLD     end do
!OLD  end do
!OLD  end do
!  write(6, "('# ORMAS: test 2')")

  p1x_beta2 = 0
  h1x_beta2 = 0
  eq1x_beta2 = 0
  sgn1x_beta2 = 0
  do istr2 = 1, nstr_beta
     istr = map2to1_beta(istr2)
     mval_beta2(istr2) = mval_beta(istr)
     dist_str_beta2(1,istr2) = dist_str_beta(1,istr)

     onv_beta2(1:nact, istr2) = onv_beta(1:nact, istr)
     orb_beta2(0:nelact(2), istr2) = orb_beta(0:nelact(2), istr)

     n1x_beta2(-3:1, istr2) = n1x_beta(-3:1, istr)
     do i1x = 1, n1x_beta2(1,istr2)
        p1x_beta2(i1x, istr2) = p1x_beta(i1x, istr)
        h1x_beta2(i1x, istr2) = h1x_beta(i1x, istr)
        sgn1x_beta2(i1x, istr2) = sgn1x_beta(i1x, istr)
        if (eq1x_beta(i1x,istr) > 0) then
           eq1x_beta2(i1x,istr2) = map1to2_beta(eq1x_beta(i1x, istr))
        else
           eq1x_beta2(i1x,istr2) = eq1x_beta(i1x, istr)
        end if
     end do
     n1x_m_beta2(-mmax2:mmax2, istr2) = n1x_m_beta(-mmax2:mmax2, istr)
     map1x_m_beta2(1:max1x_m_beta, -mmax2:mmax2, istr2) = map1x_m_beta(1:max1x_m_beta, -mmax2:mmax2, istr)
  end do
!OLD  n1xr_beta2(1:nact, 1:nact) = n1xr_beta(1:nact, 1:nact)
!OLD  r1xr_beta2 = 0
!OLD  l1xr_beta2 = 0
!OLD  sgn1xr_beta2 = 0
!OLD  do iact = 1, nact
!OLD  do jact = 1, nact
!OLD     do ii = 1, n1xr_beta2(jact,iact)
!OLD        r1xr_beta2(ii,jact,iact) = map1to2_beta(r1xr_beta(ii,jact,iact))
!OLD        l1xr_beta2(ii,jact,iact) = map1to2_beta(l1xr_beta(ii,jact,iact))
!OLD        sgn1xr_beta2(ii,jact,iact) = sgn1xr_beta(ii,jact,iact)
!OLD     end do
!OLD  end do
!OLD  end do
!  write(6, "('# ORMAS: test 3')")

  ! allowed alpha strings for a given beta strings and
  ! allowed beta strings for a given alpha strings
  llstr_alph_beta = 0
  do istr2 = 1, nstr_beta
     mstrb = mval_beta2(istr2)
     idist = dist_str_beta2(1,istr2)
     tmp = 0
     found = .false.
     do jdist = 1, ndist_alph
        if (det_allowed(jdist,idist) == 0) cycle
        do mstra = mmin_alph, mmax_alph
           if (mstra+mstrb .ne. mtot) cycle
           lls = llstr_dist_m_alph(jdist,mstra)
           uls =  nstr_dist_m_alph(jdist,mstra) + lls - 1
           do jstr2 = lls, uls
              tmp = tmp + 1
              if (.not. found) then
                 found = .true.
                 llstr_alph_beta(istr2) = jstr2
              end if
           end do
        end do
     end do
     nstr_alph_beta(istr2) = tmp
  end do
  llstr_beta_alph = 0
  do istr2 = 1, nstr_alph
     mstra = mval_alph2(istr2)
     idist = dist_str_alph2(1,istr2)
     tmp = 0
     found = .false.
     do jdist = 1, ndist_beta
        if (det_allowed(idist,jdist) == 0) cycle
        do mstrb = mmin_beta, mmax_beta
           if (mstra+mstrb .ne. mtot) cycle
           lls = llstr_dist_m_beta(jdist,mstrb)
           uls =  nstr_dist_m_beta(jdist,mstrb) + lls - 1
           do jstr2 = lls, uls
              tmp = tmp + 1
              if (.not. found) then
                 found = .true.
                 llstr_beta_alph(istr2) = jstr2
              end if
           end do
        end do
     end do
     nstr_beta_alph(istr2) = tmp
  end do

!debug
!  write(6, "('# ORMAS: test 4')")
!  do istr2 = 1, nstr_alph
!     write(6,"('strA,M,dist,llstrB(A),ulstrB(A): ', 5i10)") &
!          istr2, &
!          mval_alph2(istr2), &
!          dist_str_alph2(1,istr2), &
!          llstr_beta_alph(istr2), &
!          llstr_beta_alph(istr2)+nstr_beta_alph(istr2)-1
!  end do
!  stop
!debug

  tmp = 0
  ntot_alph_beta = 0
!##########
  do istr2 = 1, nstr_beta
     ntot_alph_beta(istr2) = tmp - llstr_alph_beta(istr2) + 1
     tmp = tmp + nstr_alph_beta(istr2)
  end do
!new  do istr2 = 1, nstr_beta
!new     ntot_alph_beta(istr2) = tmp! + 1
!new     tmp = tmp + nstr_alph_beta(istr2)
!new  end do
!##########

  tmp = 0
  ntot_beta_alph = 0
!##########
  do istr2 = 1, nstr_alph
     ntot_beta_alph(istr2) = tmp - llstr_beta_alph(istr2) + 1
     tmp = tmp + nstr_beta_alph(istr2)
  end do
!new  do istr2 = 1, nstr_alph
!new     ntot_beta_alph(istr2) = tmp! + 1
!new     tmp = tmp + nstr_beta_alph(istr2)
!new  end do
!##########
 
!debug
!  write(6, "('# ORMAS: test 5')")
!  it1 = 0
!  do istr2 = 1, nstr_beta
!     write(6,"('ntot: ',5i5,'-->',i5,':',i5,'-->',i5)") &
!          istr2, &
!          mval_beta(istr2), &
!          llstr_alph_beta(istr2), &
!          nstr_alph_beta(istr2), &
!          ntot_alph_beta(istr2)+llstr_alph_beta(istr2), &
!          ntot_alph_beta(istr2)+llstr_alph_beta(istr2)+nstr_alph_beta(istr2)-1, &
!!new          ntot_alph_beta(istr2)+1, &
!!new          ntot_alph_beta(istr2)+nstr_alph_beta(istr2), &
!          it1+1, &
!          it1+nstr_alph_beta(istr2)
!     it1 = it1 + nstr_alph_beta(istr2)
!  end do
!!  stop
!debug

  ! number of M-adapted determinants
  if (.not. cic_old) then
     ndetx = 0
     do istr2 = 1, nstr_beta
        ndetx = ndetx + nstr_alph_beta(istr2)
     end do
  else
     ndetx = nstr_alph * nstr_beta
  end if

!  allocate(mapf_detx(1:nstr_alph, 1:nstr_beta))
  allocate(mapr_detx(1:2, 1:ndetx))
!  mapf_detx = 0
  mapr_detx = 0

  ! string pair <--> determinant map
  if (.not. cic_old) then
     ndety = 0
     do istr2 = 1, nstr_beta
        lls = llstr_alph_beta(istr2)
        uls =  nstr_alph_beta(istr2) + lls - 1
        do jstr2 = lls, uls
           ndety = ndety + 1
!           mapf_detx(jstr2,istr2) = ndety
           mapr_detx(1,ndety) = jstr2
           mapr_detx(2,ndety) = istr2
        end do
     end do
     if (ndetx .ne. ndety) stop 'algorithmic failure in ormas_madapt_str (3).'
  else
     ndety = 0
     do istr2 = 1, nstr_beta
        do jstr2 = 1, nstr_alph
           ndety = ndety + 1
!           mapf_detx(jstr2,istr2) = ndety
           mapr_detx(1,ndety) = jstr2
           mapr_detx(2,ndety) = istr2
        end do
     end do
     if (ndetx .ne. ndety) stop 'algorithmic failure in ormas_madapt_str (4).'
  end if

!  !DEBUG
!  do istr2 = 1, nstr_beta
!     lls = llstr_alph_beta(istr2)
!     uls =  nstr_alph_beta(istr2) + lls - 1
!     do jstr2 = lls, uls
!        write(6, "('# ORMAS: mapf_detx ', 4i10)") jstr2,istr2,mapf_detx(jstr2,istr2),ntot_alph_beta(istr2)+jstr2
!     end do
!  end do
!  !DEBUG

  ! copy to the original variables
  if (.not. cic_old) then
     mval_alph = mval_alph2
     mval_beta = mval_beta2
     dist_str_alph = dist_str_alph2
     dist_str_beta = dist_str_beta2
     onv_alph =      onv_alph2
     onv_beta =      onv_beta2
     orb_alph =      orb_alph2
     orb_beta =      orb_beta2
     n1x_alph =      n1x_alph2
     p1x_alph =      p1x_alph2
     h1x_alph =      h1x_alph2
     eq1x_alph =     eq1x_alph2
     sgn1x_alph =    sgn1x_alph2
     n1x_m_alph =    n1x_m_alph2
     map1x_m_alph =  map1x_m_alph2
     n1x_beta =      n1x_beta2
     p1x_beta =      p1x_beta2
     h1x_beta =      h1x_beta2
     eq1x_beta =     eq1x_beta2
     sgn1x_beta =    sgn1x_beta2
     n1x_m_beta =    n1x_m_beta2
     map1x_m_beta =  map1x_m_beta2
!OLD     n1xr_alph =     n1xr_alph2
!OLD     r1xr_alph =     r1xr_alph2
!OLD     l1xr_alph =     l1xr_alph2
!OLD     sgn1xr_alph =   sgn1xr_alph2
!OLD     n1xr_beta =     n1xr_beta2
!OLD     r1xr_beta =     r1xr_beta2
!OLD     l1xr_beta =     l1xr_beta2
!OLD     sgn1xr_beta =   sgn1xr_beta2
     if (iprint > 2) then
        write(6, "('# ORMAS: alpha-strings (M-ordered)', i5)") nstr_alph
        write(6, "('      #')", advance = 'no')
        write(6, "('   mval')", advance = 'no')
        write(6, "('   map2')", advance = 'no')
        write(6, "('   dist')", advance = 'no')
        write(6, *)
        do istr = 1, nstr_alph
           write(6, "(4i7)", advance = 'no') istr, mval_alph(istr), &
                map2to1_alph(istr), dist_str_alph(1,istr)
           write(6, "(2x)", advance = 'no')
           call ormas_occvec_print(6, .false., nact, onv_alph(1,istr))
           write(6, "(2x)", advance = 'no')
           if (nelact(1) > 0) call ormas_orbvec_print(.false., nelact(1), orb_alph(1,istr))
           write(6, *)
        end do
        write(6, "('# ORMAS: beta-strings (M-ordered)', i7)") nstr_beta
        write(6, "('      #')", advance = 'no')
        write(6, "('   mval')", advance = 'no')
        write(6, "('   map2')", advance = 'no')
        write(6, "('   dist')", advance = 'no')
        write(6, *)
        do istr = 1, nstr_beta
           write(6, "(4i7)", advance = 'no') istr, mval_beta(istr), &
                map2to1_beta(istr), dist_str_beta(1,istr)
           write(6, "(2x)", advance = 'no')
           call ormas_occvec_print(6, .false., nact, onv_beta(1,istr))
           write(6, "(2x)", advance = 'no')
           if (nelact(2) > 0) call ormas_orbvec_print(.false., nelact(2), orb_beta(1,istr))
           write(6, *)
        end do
     end if
     write(6, "('# ORMAS: M-adapted determinants ', 3i20)") ndetx, ndet, nstr_alph*nstr_beta
     !debug
     if (iprint > 2) then
        do idet = 1, ndetx
           write(6, "(7i7)", advance = 'no') idet, mapr_detx(1,idet), mapr_detx(2,idet), &
                mval_alph(mapr_detx(1,idet)), mval_beta(mapr_detx(2,idet)), &
                dist_str_alph(1,mapr_detx(1,idet))-1, &
                dist_str_beta(1,mapr_detx(2,idet))-1
           write(6, "(2x)", advance = 'no')
!           if (nelact(1) > 0) call ormas_occvec_print(6,.false.,nact,onv_alph(1,mapr_detx(1,idet)))
           call ormas_occvec_print(6,.false.,nact,onv_alph(1,mapr_detx(1,idet)))
           write(6, "(' x ')", advance = 'no')
!           if (nelact(2) > 0) call ormas_occvec_print(6,.false.,nact,onv_beta(1,mapr_detx(2,idet)))
           call ormas_occvec_print(6,.false.,nact,onv_beta(1,mapr_detx(2,idet)))
           write(6, *)
        end do
     end if
     if (iprint > 2) then
        write(6, "(' # ORMAS: alpha int1x (M-ordered)')")
        call ormas_int1x_print_spin(nact, nstr_alph, onv_alph, n1x_alph, &
             & p1x_alph, h1x_alph, eq1x_alph, sgn1x_alph)
        write(6, "(' # ORMAS: beta int1x (M-ordered)')")
        call ormas_int1x_print_spin(nact, nstr_beta, onv_beta, n1x_beta, &
             & p1x_beta, h1x_beta, eq1x_beta, sgn1x_beta)
!OLD        write(6, "(' # ORMAS: alpha int1xr (M-ordered)')")
!OLD        call ormas_int1xr_print_spin(nelact(1), nact, nstr_alph, onv_alph, &
!OLD             & n1xr_alph, r1xr_alph, l1xr_alph, sgn1xr_alph)
!OLD        write(6, "(' # ORMAS: beta int1xr (M-ordered)')")
!OLD        call ormas_int1xr_print_spin(nelact(1), nact, nstr_alph, onv_alph, &
!OLD             & n1xr_alph, r1xr_alph, l1xr_alph, sgn1xr_alph)
     end if
  end if

  deallocate(mval_alph2)
  deallocate(mval_beta2)
  deallocate(dist_str_alph2)
  deallocate(dist_str_beta2)
  deallocate(onv_alph2)
  deallocate(onv_beta2)
  deallocate(orb_alph2)
  deallocate(orb_beta2)
  deallocate(n1x_alph2)
  deallocate(p1x_alph2)
  deallocate(h1x_alph2)
  deallocate(eq1x_alph2)
  deallocate(sgn1x_alph2)
  deallocate(n1x_m_alph2)
  deallocate(map1x_m_alph2)
  deallocate(n1x_beta2)
  deallocate(p1x_beta2)
  deallocate(h1x_beta2)
  deallocate(eq1x_beta2)
  deallocate(sgn1x_beta2)
  deallocate(n1x_m_beta2)
  deallocate(map1x_m_beta2)
!OLD  deallocate(n1xr_alph2)
!OLD  deallocate(r1xr_alph2)
!OLD  deallocate(l1xr_alph2)
!OLD  deallocate(sgn1xr_alph2)
!OLD  deallocate(n1xr_beta2)
!OLD  deallocate(r1xr_beta2)
!OLD  deallocate(l1xr_beta2)
!OLD  deallocate(sgn1xr_beta2)

!  write(6, "('# ORMAS: test 6')")
!  it1 = 0
!  do istr = 1, nstr_beta
!     mstrb = mval_beta(istr)
!     mstra = mtot-mstrb
!     idist = dist_str_beta(1,istr)
!     do jdist = 1, ndist_alph
!        if (det_allowed(jdist,idist) == 0) cycle
!        lla = llstr_dist_m_alph(idist,mstra)
!        ula =  nstr_dist_m_alph(idist,mstra) + lla - 1
!        it1 = it1 + nstr_dist_m_alph(idist,mstra)
!        write(6,"('ntot: ',6i5,': ',i5,'-->',i5)") &
!             istr,mstrb,mstra,idist,jdist,it1, &
!             ntot_alph_beta(istr)+lla, &
!             ntot_alph_beta(istr)+ula
!     end do
!  end do
!  stop "STOP for debug @ ormas_madapt_str."

end subroutine ormas_madapt_str
!################################################################################
!################################################################################
