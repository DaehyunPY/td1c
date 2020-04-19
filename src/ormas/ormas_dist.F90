!################################################################################
subroutine ormas_dist()

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : iprint
  use mod_ormas, only : nelact, nsub, min_sub, max_sub
  use mod_ormas, only : min_sub_alph, max_sub_alph, min_sub_beta, max_sub_beta
  use mod_ormas, only : ndist, ndist_alph, ndist_beta
  use mod_ormas, only : dist, dist_alph, dist_beta

  implicit none
  integer(c_int) :: isub, idist

!old  call ormas_dist_itr(.false., nelact(1), nsub, min_sub_alph, max_sub_alph, ndist_alph)
!old  call ormas_dist_itr(.false., nelact(2), nsub, min_sub_beta, max_sub_beta, ndist_beta)
!old  call ormas_dist_itr(.false., nelact(3), nsub, min_sub, max_sub, ndist)
  call ormas_dist_ndist(nelact(1), nsub, min_sub_alph, max_sub_alph, ndist_alph)
  call ormas_dist_ndist(nelact(2), nsub, min_sub_beta, max_sub_beta, ndist_beta)
  call ormas_dist_ndist(nelact(3), nsub, min_sub, max_sub, ndist)

  allocate(dist_alph(1:nsub, 1:ndist_alph))
  allocate(dist_beta(1:nsub, 1:ndist_beta))
  allocate(dist(1:nsub, 1:ndist))
  dist_alph(1:nsub, 1:ndist_alph) = 0
  dist_beta(1:nsub, 1:ndist_beta) = 0
  dist(1:nsub, 1:ndist) = 0

  call ormas_dist_itr(nelact(1), nsub, min_sub_alph, max_sub_alph, &
       & ndist_alph, dist_alph)
  call ormas_dist_itr(nelact(2), nsub, min_sub_beta, max_sub_beta, &
       & ndist_beta, dist_beta)
  call ormas_dist_itr(nelact(3), nsub, min_sub, max_sub, ndist, dist)

  if (iprint > 0) then
     write(6, "('# ORMAS: alpha-distribution', i5)") ndist_alph
     write(6, "(10x)", advance = 'no')
     do isub = 1, nsub
        write(6, "(i10)", advance = 'no') isub
     end do
     write(6, *)
     do idist = 1, ndist_alph
        write(6, "(i10)", advance = 'no') idist
        do isub = 1, nsub
           write(6, "(i10)", advance = 'no') dist_alph(isub, idist)
        end do
        write(6, *)
     end do

     write(6, "('# ORMAS: beta-distribution', i5)") ndist_beta
     write(6, "(10x)", advance = 'no')
     do isub = 1, nsub
        write(6, "(i10)", advance = 'no') isub
     end do
     write(6, *)
     do idist = 1, ndist_beta
        write(6, "(i10)", advance = 'no') idist
        do isub = 1, nsub
           write(6, "(i10)", advance = 'no') dist_beta(isub, idist)
        end do
        write(6, *)
     end do
  end if

end subroutine ormas_dist
!################################################################################
!################################################################################
subroutine ormas_dist_ndist(nel, nsub, min_sub, max_sub, ndist)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : iprint

  implicit none
  integer(c_int), intent(in) :: nel, nsub, min_sub(1:*), max_sub(1:*)
  integer(c_int), intent(out) :: ndist

  logical :: getnew
  integer(c_int) :: isub, nel_dyn, rest
  integer(c_int), allocatable :: max_dyn(:)
  integer(c_int), allocatable :: dist1(:)
  integer(c_int), allocatable :: dist1_dyn(:)
  logical(c_bool), external :: ormas_dist_chk_dplus

  allocate(max_dyn(1:nsub))
  allocate(dist1(1:nsub))
  allocate(dist1_dyn(1:nsub))

  ! dynamical boundaries
  nel_dyn = nel
  do isub = 1, nsub
     nel_dyn = nel_dyn - min_sub(isub)
     max_dyn(isub) = max_sub(isub) - min_sub(isub)
  end do

  ! first distribution
  rest = nel_dyn
  dist1_dyn(1) = min(max_dyn(1), rest)
  do isub = 2, nsub
     rest = rest - dist1_dyn(isub - 1)
     dist1_dyn(isub) = min(max_dyn(isub), rest)
  end do

  ! recursive search of possible distributions
  ndist = 0
  getnew = .true.
  do while(getnew)
     dist1(1:nsub) = dist1_dyn(1:nsub) + min_sub(1:nsub)
     if (ormas_dist_chk_dplus(nel,min_sub,max_sub,dist1)) ndist = ndist + 1   
     call ormas_dist_next(nel_dyn, nsub, max_dyn, getnew, dist1_dyn)
  end do

  deallocate(dist1_dyn)
  deallocate(dist1)
  deallocate(max_dyn)

end subroutine ormas_dist_ndist
!################################################################################
!################################################################################
subroutine ormas_dist_itr(nel, nsub, min_sub, max_sub, ndist, dist)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : iprint

  implicit none
  integer(c_int), intent(in) :: nel, nsub, min_sub(1:*), max_sub(1:*)
  integer(c_int), intent(out) :: ndist
  integer(c_int), intent(out) :: dist(1:nsub, 1:*)

  logical :: getnew
  integer(c_int) :: isub, nel_dyn, rest
  integer(c_int), allocatable :: max_dyn(:)
  integer(c_int), allocatable :: dist1(:)
  integer(c_int), allocatable :: dist1_dyn(:)
  logical(c_bool), external :: ormas_dist_chk_dplus

  allocate(max_dyn(1:nsub))
  allocate(dist1(1:nsub))
  allocate(dist1_dyn(1:nsub))

  ! dynamical boundaries
  nel_dyn = nel
  do isub = 1, nsub
     nel_dyn = nel_dyn - min_sub(isub)
     max_dyn(isub) = max_sub(isub) - min_sub(isub)
  end do

  ! first distribution
  rest = nel_dyn
  dist1_dyn(1) = min(max_dyn(1), rest)
  do isub = 2, nsub
     rest = rest - dist1_dyn(isub - 1)
     dist1_dyn(isub) = min(max_dyn(isub), rest)
  end do

  ! recursive search of possible distributions
  ndist = 0
  getnew = .true.
  do while(getnew)
     dist1(1:nsub) = dist1_dyn(1:nsub) + min_sub(1:nsub)
     if (ormas_dist_chk_dplus(nel,min_sub,max_sub,dist1)) then
        ndist = ndist + 1
        dist(1:nsub, ndist) = dist1(1:nsub)
     end if
     call ormas_dist_next(nel_dyn, nsub, max_dyn, getnew, dist1_dyn)
  end do

  deallocate(dist1_dyn)
  deallocate(dist1)
  deallocate(max_dyn)

end subroutine ormas_dist_itr
!################################################################################
!################################################################################
subroutine ormas_dist_next(nel, nsub, max_sub, getnew, dist)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : iprint

  implicit none
  integer(c_int), intent(in) :: nel, nsub, max_sub(1:nsub)
  logical, intent(inout) :: getnew
  integer(c_int), intent(inout) :: dist(1:nsub)

  integer(c_int) :: isub, jsub, ksub, rest

  getnew = .false.
  do isub = nsub, 1, -1
     if (dist(isub) > 0) then
        do jsub = isub + 1, nsub
           if (dist(jsub) < max_sub(jsub)) then
              dist(isub) = dist(isub) - 1
              rest = nel - sum(dist(1:(jsub-1)))
              do ksub = jsub, nsub
                 dist(ksub) = min(max_sub(ksub), rest)
                 rest = rest - dist(ksub)
!bug                 if (rest == 0) exit
              end do
              getnew = .true.
              return
           end if
        end do
     end if
  end do

end subroutine ormas_dist_next
!################################################################################
logical(c_bool) function ormas_dist_chk_dplus(nel, min_sub, max_sub, dist)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : iprint,dplus,dplus_type

  implicit none
  integer(c_int), intent(in) :: nel, min_sub(1:*), max_sub(1:*)
  integer(c_int), intent(in) :: dist(1:*)

  logical(c_bool) :: check
!old  integer(c_int) :: dist1, dist2, dist3, dist4

  check = .true.
  if (.not.dplus .or. dplus_type < 0) then
     continue
  else if (dplus_type == 0) then
     if ((dist(1)+dist(2).lt.nel-2 .and. dist(1).ne.max_sub(1)) .or. &
         (dist(3)+dist(4).gt.2 .and. dist(4).ne.0)) check = .false.
  else if (dplus_type == 1) then
     if (dist(2)+dist(3).gt.2 .and. dist(3).ne.0) check = .false.
  else if (dplus_type == 2) then
     if (dist(1)+dist(2).lt.nel-2 .and. dist(1).ne.max_sub(1)) check = .false.
  end if
!old  if (.not.dplus) then
!old     continue
!old  else if (dplus_type == 0) then
!old     dist1 = dist(1) + min_sub(1)
!old     dist2 = dist(2) + min_sub(2)
!old     dist3 = dist(3) + min_sub(3)
!old     dist4 = dist(4) + min_sub(4)
!old     if ((dist1+dist2.lt.nel-2 .and. dist1.ne.max_sub(1)) .or. &
!old         (dist3+dist4.gt.2 .and. dist4.ne.0)) check = .false.
!old  else if (dplus_type == 1) then
!old     dist2 = dist(2) + min_sub(2)
!old     dist3 = dist(3) + min_sub(3)
!old     if (dist2+dist3.gt.2 .and. dist3.ne.0) check = .false.
!old  else if (dplus_type == 2) then
!old     dist1 = dist(1) + min_sub(1)
!old     dist2 = dist(2) + min_sub(2)
!old     if (dist1+dist2.lt.nel-2 .and. dist1.ne.max_sub(1)) check = .false.
!old  end if

  ormas_dist_chk_dplus = check
  return

end function ormas_dist_chk_dplus
!################################################################################
integer(c_int) function ormas_dist_type_dplus(nel, min_sub, max_sub, dist)
!
! type = -1 for Dplus violating term
! type = +1 for |0,-3|3,0>
! type = +2 for |-3|3,0>
! type = +3 for |0,-3|3>
! type =  0 otherwise
!
  use, intrinsic :: iso_c_binding
  use mod_ormas, only : iprint,dplus,dplus_type

  implicit none
  integer(c_int), intent(in) :: nel, min_sub(1:*), max_sub(1:*)
  integer(c_int), intent(in) :: dist(1:*)

  integer(c_int) :: type
!old  integer(c_int) :: dist2, dist3
  logical(c_bool), external :: ormas_dist_chk_dplus

  type = 0
  if (.not.dplus .or. dplus_type < 0) then
     continue
  else if (.not.ormas_dist_chk_dplus(nel,min_sub,max_sub,dist)) then
     type = -1
  else if (dplus_type == 0) then
     if (dist(3).gt.2) type = 1
  else if (dplus_type == 1) then
     if (dist(2).gt.2) type = 1
  else if (dplus_type == 2) then
     if (dist(3).gt.2) type = 1
  end if
!old  if (.not.dplus) then
!old     continue
!old  else if (.not.ormas_dist_chk_dplus(nel,min_sub,max_sub,dist)) then
!old     type = -1
!old  else if (dplus_type == 0) then
!old     dist3 = dist(3) + min_sub(3)
!old     if (dist3.gt.2) type = 1
!old  else if (dplus_type == 1) then
!old     dist2 = dist(2) + min_sub(2)
!old     if (dist2.gt.2) type = 1
!old  else if (dplus_type == 2) then
!old     dist3 = dist(3) + min_sub(3)
!old     if (dist3.gt.2) type = 1
!old  end if

  ormas_dist_type_dplus = type
  return

end function ormas_dist_type_dplus
!################################################################################
