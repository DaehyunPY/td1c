!################################################################################
integer(c_int) function ormas_allowed(idist_alph, idist_beta)

  use, intrinsic :: iso_c_binding

  use mod_ormas, only : nsub
  use mod_ormas, only : ndist
  use mod_ormas, only : dist, dist_alph, dist_beta

  implicit none
  integer(c_int), intent(in) :: idist_alph, idist_beta

  logical :: found
  integer(c_int) :: isub, idist, neld, nelx

  ormas_allowed = -1

  do idist = 1, ndist
     found = .true.
     do isub = 1, nsub
        neld = dist(isub, idist)
        nelx = dist_alph(isub, idist_alph) + dist_beta(isub, idist_beta)
        if (nelx /= neld) found = .false.
     end do
     if (found) then
        ormas_allowed = idist
        exit
     end if
  end do

! 2019/6/13, commented out
  if (.not. found) stop 'error in ormas_allowed.'

  return

end function ormas_allowed
!################################################################################
