!################################################################################
integer(c_int) function ormas_elec2fun(iel,istr,nela,orb)
  use,intrinsic :: iso_c_binding
  use mod_ormas,only : ncore

  implicit none
  !--------------------------------------------------------------------
  integer(c_int),intent(in) :: iel,istr,nela,orb(0:nela,*)
  !--------------------------------------------------------------------

  if (iel <= ncore) then
     ormas_elec2fun = iel
  else
     ormas_elec2fun = ncore+orb(iel-ncore,istr)
  end if
  return

end function ormas_elec2fun
!################################################################################
