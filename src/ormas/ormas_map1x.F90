!#######################################################################
integer(c_int) function ormas_map1x(nelec, norb, ndist, dist, nstr_dist, &
     & nstr_dist_sub, arc, occvec)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nact, nsub, norb_sub, lorb_sub

  implicit none
  integer(c_int), intent(in) :: nelec, norb, ndist
  integer(c_int), intent(in) :: dist(1:nsub, 1:ndist)
  integer(c_int), intent(in) :: nstr_dist(1:ndist)
  integer(c_int), intent(in) :: nstr_dist_sub(1:nsub, 1:ndist)
  integer(c_int), intent(in) :: arc(1:nact, 1:*)
  integer(c_int), intent(in) :: occvec(1:*)

  integer(c_int), allocatable :: disvec(:)
  integer(c_int), allocatable :: substr(:)
  logical :: found
  integer(c_int) :: xdist, offset, nel, llorb, ulorb, rest
  integer(c_int) :: idist, isub, jsub, iorb, imap
  integer(c_int), external :: ormas_map1x_sub

  if (nelec == 0) then
     ormas_map1x = 1
     return
  end if

  allocate(disvec(1:nsub))
  allocate(substr(1:nsub))

  ! subspace occupations
  disvec(1:nsub) = 0
  do isub = 1, nsub
     do iorb = lorb_sub(1, isub), lorb_sub(2, isub)
        disvec(isub) = disvec(isub) + min(1, occvec(iorb))
     end do
  end do

  ! which distributin?
  xdist = -1
  offset = 0
  do idist = 1, ndist
     found = .true.
     do isub = 1, nsub
        found = found .and. (disvec(isub) == dist(isub, idist))
     end do
     if (found) then
        xdist = idist
        exit
     else
        offset = offset + nstr_dist(idist)
     end if
  end do

  if (xdist < 0) stop 'error in ormas_map1x.'

  ! short-code numberings
  do isub = 1, nsub
     nel = disvec(isub)
     llorb = lorb_sub(1, isub)
     ulorb = lorb_sub(2, isub)
     substr(isub) = ormas_map1x_sub(nel, norb_sub(isub), arc, occvec(llorb:ulorb))
  end do

  ! global numvering
  imap = offset + 1
  do isub = 1, nsub
     rest = 1
     do jsub = isub + 1, nsub
        rest = rest * nstr_dist_sub(jsub, xdist)
     end do
     imap = imap + (substr(isub) - 1) * rest
  end do

  ormas_map1x = imap

  deallocate(substr)
  deallocate(disvec)

end function ormas_map1x
!#######################################################################
!#######################################################################
integer(c_int) function ormas_map1x_sub(nelec, norb, arc, occvec)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nact

  implicit none
  integer(c_int), intent(in) :: nelec, norb
  integer(c_int), intent(in) :: arc(1:nact, 1:*)
  integer(c_int), intent(in) :: occvec(1:*)

  integer(c_int), allocatable :: string(:)
  integer(c_int) :: iorb, iel, nel, imap

  if (nelec == 0) then
     ormas_map1x_sub = 1
     return
  end if

  allocate(string(1:norb))

  nel = 0
  do iorb = 1, norb
     if (occvec(iorb) /= 0) then
        nel = nel + 1
        string(nel) = iorb
     end if
  end do

  imap = 1
  do iel = 1, nel
     imap = imap + arc(string(iel),iel)
  end do
  ormas_map1x_sub = imap

  deallocate(string)

end function ormas_map1x_sub
!#######################################################################
