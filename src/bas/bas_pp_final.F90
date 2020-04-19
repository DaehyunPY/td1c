!######################################################################
subroutine bas_pp_final()

  use, intrinsic :: iso_c_binding
  use mod_control, only : PSP
  use mod_bas, only : pp_vloc, pp_vlocHF, pp_pproj, pp_fproj, pp_gproj

  implicit none
  if (.not. PSP) return
  deallocate(pp_vloc)
  deallocate(pp_vlocHF)
  deallocate(pp_pproj)
  deallocate(pp_fproj)
  deallocate(pp_gproj)

end subroutine bas_pp_final
!######################################################################
