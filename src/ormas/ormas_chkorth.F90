!################################################################################
logical(c_bool) function ormas_chkorth(i1,j1,istr,jstr,nela,orb)
  use,intrinsic :: iso_c_binding
  use mod_ormas,only : ncore

  implicit none
  !--------------------------------------------------------------------
  integer(c_int),intent(in) :: i1,j1,istr,jstr,nela,orb(0:nela,*)
  !--------------------------------------------------------------------

  ormas_chkorth = &
       i1<=ncore.and.i1==j1 .or. &
       i1>ncore.and.orb(i1-ncore,istr)==orb(j1-ncore,jstr)
  return

end function ormas_chkorth
!################################################################################
