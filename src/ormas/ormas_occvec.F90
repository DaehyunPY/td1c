!################################################################################
subroutine ormas_occvec(nel, norb, istr, arc, occvec)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : iprint
  use mod_ormas, only : nact

  implicit none
  integer(c_long), intent(in) :: nel, norb, istr
  integer(c_long), intent(in) :: arc(1:nact, 1:*)
  integer(c_long), intent(out) :: occvec(1:*)
  integer(c_long) :: iel, iorb, rest

  occvec(1:norb) = 0

  iel = nel
  rest = istr - 1
  do iorb = norb, 1, -1
     if (iel <= 0) exit
     if (rest >= arc(iorb, iel)) then
        occvec(iorb) = 1
        rest = rest - arc(iorb, iel)
        iel = iel - 1
     end if
  end do

!debug
!  do iorb = 1, norb
!     write(6, "(i1)", advance = 'no') occvec(iorb)
!  end do
!  write(6, *)
!debug

end subroutine ormas_occvec
!################################################################################
!#######################################################################
subroutine ormas_occvec_print(iow, adv, nact, occvec)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nsub, lorb_sub

  implicit none
  integer(c_long), intent(in) :: iow
  logical, intent(in) :: adv
  integer(c_long), intent(in) :: nact
  integer(c_long), intent(in) :: occvec(1:nact)

  integer(c_long) :: isub, iorb, occn

  call ormas_core_print(iow, .false.)
!  write(iow, "('|')", advance='no')

  do isub = 1, nsub - 1
     do iorb = lorb_sub(1, isub), lorb_sub(2, isub)
        occn = min(1,abs(occvec(iorb)))
        write(iow,"(i1)",advance='no') occn
     end do
     write(iow, "('.')", advance='no')
  end do

  do iorb = lorb_sub(1, nsub), lorb_sub(2, nsub)
     occn = min(1,abs(occvec(iorb)))
     write(iow,"(i1)",advance='no') occn
  end do

  if (adv) write(iow,*)

end subroutine ormas_occvec_print
!#######################################################################
!#######################################################################
subroutine ormas_core_print(iow, adv)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nfcore, ndcore

  implicit none
  integer(c_long), intent(in) :: iow
  logical, intent(in) :: adv
  integer(c_long) :: ifun

  do ifun = 1, nfcore
     write(iow, "(a1)", advance='no') 'F'
  end do
  do ifun = 1, ndcore
     write(iow, "(a1)", advance='no') 'C'
  end do
  if (adv) write(iow, *)

end subroutine ormas_core_print
!#######################################################################
