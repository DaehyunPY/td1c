!################################################################################
subroutine ormas_orbvec(norb, onv, orb)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : iprint

  implicit none
  integer(c_long), intent(in) :: norb
  integer(c_long), intent(in) :: onv(1:*)
  integer(c_long), intent(out) :: orb(0:*)

  integer(c_long) :: nel_tot, iorb

  nel_tot = 0
  do iorb = 1, norb
     if (onv(iorb) /= 0) then
        nel_tot = nel_tot + 1
        orb(nel_tot) = iorb
     end if
  end do

  if (nel_tot > 0) then
     orb(0) = 1
  else
     orb(0) = 0
  end if

end subroutine ormas_orbvec
!################################################################################
!#######################################################################
subroutine ormas_orbvec_print(adv, nel, orbvec)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : ncore

  implicit none
  logical, intent(in) :: adv
  integer(c_long), intent(in) :: nel
  integer(c_long), intent(in) :: orbvec(1:nel)
  integer(c_long) :: iel

  do iel = 1, ncore - 1
     write(6,"(i2,',')",advance='no') iel
  end do
  if (ncore > 0) write(6,"(i2)",advance='no') ncore

  if (nel > 0) then
     if (ncore > 0) write(6,"(',')",advance='no')
     do iel = 1, nel - 1
        write(6,"(i2,',')",advance='no') ncore + orbvec(iel)
     end do
     write(6,"(i2)",advance='no') ncore + orbvec(nel)
  end if

  if (adv) write(6,*)

end subroutine ormas_orbvec_print
!#######################################################################
