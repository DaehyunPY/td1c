!################################################################################
subroutine ormas_nstr

  use, intrinsic :: iso_c_binding

  use mod_ormas, only : iprint
  use mod_ormas, only : nelact, nsub, norb_sub, min_sub, max_sub, ncore
  use mod_ormas, only : min_sub_alph, max_sub_alph, min_sub_beta, max_sub_beta
  use mod_ormas, only : ndist, ndist_alph, ndist_beta, ormas_donly, ormas_sd1, read_allowed
  use mod_ormas, only : dist, dist_alph, dist_beta, ndet_dist
  use mod_ormas, only : ndet, det_allowed, dplus
  use mod_ormas, only : nstr_alph, nstr_alph_dist, lstr_alph_dist, nstr_alph_dist_sub
  use mod_ormas, only : nstr_beta, nstr_beta_dist, lstr_beta_dist, nstr_beta_dist_sub

  implicit none
  logical(c_bool) :: allowed
  integer(c_int) :: idtot, isub, idist, jdist, istr, jstr, idet, nstrx, sumab, refab, &
       difa, difb, difab, tmp
  integer(c_int), allocatable :: dist1(:)
  integer(c_int), external :: util_bicoeff
  integer(c_int), external :: ormas_allowed
  logical(c_bool), external :: ormas_dist_chk_dplus

  allocate(det_allowed(1:ndist_alph, 1:ndist_beta))
  allocate(nstr_alph_dist_sub(1:nsub, 1:ndist_alph))
  allocate(nstr_beta_dist_sub(1:nsub, 1:ndist_beta))
  allocate(nstr_alph_dist(1:ndist_alph))
  allocate(nstr_beta_dist(1:ndist_beta))
  allocate(lstr_alph_dist(1:2, 1:ndist_alph))
  allocate(lstr_beta_dist(1:2, 1:ndist_beta))

  ! number of strings
  tmp = 0
  do idist = 1, ndist_alph
     nstr_alph_dist(idist) = 1
     do isub = 1, nsub
        nstrx = util_bicoeff(norb_sub(isub), dist_alph(isub, idist))
        nstr_alph_dist(idist) = nstr_alph_dist(idist) * nstrx
        nstr_alph_dist_sub(isub, idist) = nstrx
     end do
     lstr_alph_dist(1, idist) = tmp + 1
     lstr_alph_dist(2, idist) = tmp + nstr_alph_dist(idist)
     tmp = lstr_alph_dist(2, idist)
  end do
  nstr_alph = sum(nstr_alph_dist(1:ndist_alph))

  tmp = 0
  do idist = 1, ndist_beta
     nstr_beta_dist(idist) = 1
     do isub = 1, nsub
        nstrx = util_bicoeff(norb_sub(isub), dist_beta(isub, idist))
        nstr_beta_dist(idist) = nstr_beta_dist(idist) * nstrx
        nstr_beta_dist_sub(isub, idist) = nstrx
     end do
     lstr_beta_dist(1, idist) = tmp + 1
     lstr_beta_dist(2, idist) = tmp + nstr_beta_dist(idist)
     tmp = lstr_beta_dist(2, idist)
  end do
  nstr_beta = sum(nstr_beta_dist(1:ndist_beta))

  ! allowed combinations of alpha and beta strings
  if (read_allowed) then
!    call ormas_allowed_read
     stop "ormas_allowed_read nyi."
  else
     det_allowed(1:ndist_alph, 1:ndist_beta) = 0
  end if

  allocate(dist1(1:nsub))
  do idist = 1, ndist_beta
     do jdist = 1, ndist_alph
        if (det_allowed(jdist, idist) < 0) then
           det_allowed(jdist, idist) = 0
           cycle
        end if

        allowed = .true.
        do isub = 1, nsub
           sumab = dist_alph(isub, jdist) + dist_beta(isub, idist)
           if ((min_sub(isub) > sumab) .or. &
             & (max_sub(isub) < sumab)) then
              allowed = .false.
              exit
           end if
        end do

! 20180208        !##### DONLY #####
! 20180208        if (ormas_donly) then
! 20180208           do isub = 1, nsub
! 20180208              refab = dist_alph(isub, 1) + dist_beta(isub, 1)
! 20180208              sumab = dist_alph(isub, jdist) + dist_beta(isub, idist)
! 20180208              difab = abs(sumab - refab)
! 20180208              if (difab /= 0 .and. difab /= 2) then
! 20180208                 allowed = .false.
! 20180208                 exit
! 20180208              end if
! 20180208           end do
! 20180208        end if
! 20180208        !##### DONLY #####

        !##### S+D(SS) #####
        if (ormas_sd1) then
           do isub = 1, nsub
              difa = abs(dist_alph(isub, jdist) - dist_alph(isub, 1))
              difb = abs(dist_beta(isub, idist) - dist_beta(isub, 1))
              if (difa > 1 .or. difb > 1) then
                 allowed = .false.
                 exit
              end if
           end do
        end if

        !##### Dplus #####
        if (dplus) then
           do isub = 1, nsub
              dist1(isub) = dist_alph(isub, jdist) + dist_beta(isub, idist)
           end do
           allowed = allowed.and.ormas_dist_chk_dplus(nelact(3),min_sub,max_sub,dist1)
        end if

        !if (allowed) det_allowed(jdist, idist) = 1
        if (allowed) det_allowed(jdist, idist) = ormas_allowed(jdist, idist)
     end do
  end do
  deallocate(dist1)

  ! number of determinants
  ndet = 0
  allocate(ndet_dist(1:ndist))
  ndet_dist(1:ndist) = 0
  do idist = 1, ndist_beta
     do jdist = 1, ndist_alph
        !2019/6/13
        !if (det_allowed(jdist, idist) /= 0) then
        if (det_allowed(jdist, idist) > 0) then
           tmp = nstr_alph_dist(jdist) * nstr_beta_dist(idist)
           ndet = ndet + tmp
           idtot = det_allowed(jdist, idist)
           ndet_dist(idtot) = ndet_dist(idtot) + tmp
        end if
     end do
  end do

  if (iprint > 0) then
     write(6, "('# ORMAS: distribution', i5)") ndist
     write(6, "(10x)", advance = 'no')
     do isub = 1, nsub
        write(6, "(i10)", advance = 'no') isub
     end do
     write(6, *)
     do idist = 1, ndist
        write(6, "(i10)", advance = 'no') idist
        do isub = 1, nsub
           write(6, "(i10)", advance = 'no') dist(isub, idist)
        end do
        write(6, "(i10)") ndet_dist(idist)
     end do
     write(6, "('# total:  ')", advance = 'no')
     do isub = 1, nsub
        write(6, "(10x)", advance = 'no')
     end do
     write(6, "(i10)") ndet
  end if

  deallocate(ndet_dist)

!old  ! string pair <--> determinant map
!old  allocate(mapf_det(1:nstr_alph, 1:nstr_beta))
!old  allocate(mapr_det(1:2, 1:ndet))
!old  idet = 0
!old  mapf_det(1:nstr_alph, 1:nstr_beta) = -1
!old  do idist = 1, ndist_beta
!old     do jdist = 1, ndist_alph
!old        if (det_allowed(jdist, idist) == 0) cycle
!old        do istr = lstr_beta_dist(1, idist), lstr_beta_dist(2, idist)
!old           do jstr = lstr_alph_dist(1, jdist), lstr_alph_dist(2, jdist)
!old              idet = idet + 1
!old              mapf_det(jstr, istr) = idet
!old              mapr_det(1, idet) = jstr
!old              mapr_det(2, idet) = istr
!old           end do
!old        end do
!old     end do
!old  end do

!  write(6, "('# ORMAS: number of determinant', i10)") ndet
  if (iprint > 0) then
     write(6, "('# ORMAS: number of alpha-strings', i5)") nstr_alph
     write(6, "(40x)", advance = 'no')
     do isub = 1, nsub
        write(6, "(i10)", advance = 'no') isub
     end do
     write(6, *)
     do idist = 1, ndist_alph
        write(6, "(4i10)", advance = 'no') idist, &
             & nstr_alph_dist(idist), lstr_alph_dist(1:2, idist)
        do isub = 1, nsub
           write(6, "(i10)", advance = 'no') nstr_alph_dist_sub(isub, idist)
        end do
        write(6, *)
     end do

     write(6, "('# ORMAS: number of beta-strings', i5)") nstr_beta
     write(6, "(40x)", advance = 'no')
     do isub = 1, nsub
        write(6, "(i10)", advance = 'no') isub
     end do
     write(6, *)
     do idist = 1, ndist_beta
        write(6, "(4i10)", advance = 'no') idist, &
             & nstr_beta_dist(idist), lstr_beta_dist(1:2, idist)
        do isub = 1, nsub
           write(6, "(i10)", advance = 'no') nstr_beta_dist_sub(isub, idist)
        end do
        write(6, *)
     end do

     write(6, "('# ORMAS: det_allowed')")
     write(6, "(5x)", advance = 'no')
     do idist = 1, ndist_beta
        write(6, "(i5)", advance = 'no') idist
     end do
     write(6, *)
     do jdist = 1, ndist_alph
        write(6, "(i5)", advance = 'no') jdist
        do idist = 1, ndist_beta
           write(6, "(i5)", advance = 'no') det_allowed(jdist, idist)
        end do
        write(6, *)
     end do
  end if

end subroutine ormas_nstr
!################################################################################
