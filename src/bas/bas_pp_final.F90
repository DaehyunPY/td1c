!######################################################################
subroutine bas_pp_final()

  use, intrinsic :: iso_c_binding
  use mod_control, only : PSP, psp_type
  use mod_bas, only : pp_vloc, pp_vlocHF, pp_pproj, pp_fproj, pp_gproj
  use mod_ppatom
  use mod_oncvpsp
  use mod_tmpsp

  implicit none
  if (.not. PSP) return
  deallocate(pp_vloc)
  deallocate(pp_vlocHF)
  deallocate(pp_pproj)
  deallocate(pp_fproj)
  deallocate(pp_gproj)

  if (psp_type == 1) then
!nyi  else if (psp_type == 2) then
!nyi  else if (psp_type == 3) then
!nyi  else if (psp_type == 4) then
!nyi     call ppatom_final
!nyi  else if (psp_type == 5) then
!nyi     call ppatom_final
  else if (psp_type == 6) then
     call oncvpsp_final
  else if (psp_type == 7) then
     call tmpsp_final
  else
     stop 'bas_pp_init: bad psp_type.'
  end if

end subroutine bas_pp_final
!######################################################################
