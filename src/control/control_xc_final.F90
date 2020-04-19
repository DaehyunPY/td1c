!######################################################################
subroutine control_xc_final()

  ! eXchange and correlation

  use, intrinsic :: iso_c_binding
  use xc_f03_lib_m
  use mod_control, only : xc_funcx, xc_funcc

  implicit none

  call xc_f03_func_end(xc_funcx)
  call xc_f03_func_end(xc_funcc)

end subroutine control_xc_final
!######################################################################
