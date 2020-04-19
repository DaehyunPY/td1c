!################################################################################
subroutine ormas_dist()

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : iprint
  use mod_ormas, only : nelact, nsub, min_sub, max_sub
  use mod_ormas, only : min_sub_alph, max_sub_alph, min_sub_beta, max_sub_beta
  use mod_ormas, only : ndist, ndist_alph, ndist_beta
  use mod_ormas, only : dist, dist_alph, dist_beta

  implicit none
  integer(c_long) :: isub, idist

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
!     write(6, "('# ORMAS: distribution', i5)") ndist
!     write(6, "(10x)", advance = 'no')
!     do isub = 1, nsub
!        write(6, "(i10)", advance = 'no') isub
!     end do
!     write(6, *)
!     do idist = 1, ndist
!        write(6, "(i10)", advance = 'no') idist
!        do isub = 1, nsub
!           write(6, "(i10)", advance = 'no') dist(isub, idist)
!        end do
!        write(6, *)
!     end do

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
  integer(c_long), intent(in) :: nel, nsub, min_sub(1:*), max_sub(1:*)
  integer(c_long), intent(out) :: ndist

  logical :: getnew
  integer(c_long) :: isub, nel_dyn, rest
  integer(c_long), allocatable :: max_dyn(:)
  integer(c_long), allocatable :: dist1(:)

  allocate(max_dyn(1:nsub))
  allocate(dist1(1:nsub))

  ! dynamical boundaries
  nel_dyn = nel
  do isub = 1, nsub
     nel_dyn = nel_dyn - min_sub(isub)
     max_dyn(isub) = max_sub(isub) - min_sub(isub)
  end do

  ! first distribution
  rest = nel_dyn
  dist1(1) = min(max_dyn(1), rest)
  do isub = 2, nsub
     rest = rest - dist1(isub - 1)
     dist1(isub) = min(max_dyn(isub), rest)
  end do

  ! recursive search of possible distributions
  ndist = 0
  getnew = .true.
  do while(getnew)
     ndist = ndist + 1
     call ormas_dist_next(nel_dyn, nsub, max_dyn, getnew, dist1)
  end do

  deallocate(max_dyn)
  deallocate(dist1)

end subroutine ormas_dist_ndist
!################################################################################
!################################################################################
subroutine ormas_dist_itr(nel, nsub, min_sub, max_sub, ndist, dist)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : iprint

  implicit none
  integer(c_long), intent(in) :: nel, nsub, min_sub(1:*), max_sub(1:*)
  integer(c_long), intent(out) :: ndist
  integer(c_long), intent(out) :: dist(1:nsub, 1:*)

  logical :: getnew
  integer(c_long) :: isub, nel_dyn, rest
  integer(c_long), allocatable :: max_dyn(:)
  integer(c_long), allocatable :: dist1(:)

  allocate(max_dyn(1:nsub))
  allocate(dist1(1:nsub))

  ! dynamical boundaries
  nel_dyn = nel
  do isub = 1, nsub
     nel_dyn = nel_dyn - min_sub(isub)
     max_dyn(isub) = max_sub(isub) - min_sub(isub)
  end do

  ! first distribution
  rest = nel_dyn
  dist1(1) = min(max_dyn(1), rest)
  do isub = 2, nsub
     rest = rest - dist1(isub - 1)
     dist1(isub) = min(max_dyn(isub), rest)
  end do

  ! recursive search of possible distributions
  ndist = 0
  getnew = .true.
  do while(getnew)
     ndist = ndist + 1
     dist(1:nsub, ndist) = dist1(1:nsub) + min_sub(1:nsub)
     call ormas_dist_next(nel_dyn, nsub, max_dyn, getnew, dist1)
  end do

  deallocate(max_dyn)
  deallocate(dist1)

end subroutine ormas_dist_itr
!################################################################################
!################################################################################
subroutine ormas_dist_next(nel, nsub, max_sub, getnew, dist)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : iprint

  implicit none
  integer(c_long), intent(in) :: nel, nsub, max_sub(1:nsub)
  logical, intent(inout) :: getnew
  integer(c_long), intent(inout) :: dist(1:nsub)

  integer(c_long) :: isub, jsub, ksub, rest

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
