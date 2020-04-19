!######################################################################
subroutine control_xc_init(dft_type)

  ! eXchange and correlation

  use, intrinsic :: iso_c_binding
  use xc_f03_lib_m
  use mod_control, only : xc_funcx, xc_funcc, xc_infox, xc_infoc, func_idx, func_idc

  implicit none
  integer(c_long), intent(in) :: dft_type
  integer(c_int) :: vmajor, vminor, vmicro
  character(len=120) :: xname, cname


  func_idx = mod(dft_type, 10)
  func_idc = dft_type / 10

  call xc_f03_func_init(xc_funcx, func_idx, XC_UNPOLARIZED)
  xc_infox = xc_f03_func_get_info(xc_funcx)

  call xc_f03_func_init(xc_funcc, func_idc, XC_UNPOLARIZED)
  xc_infoc = xc_f03_func_get_info(xc_funcc)

  call xc_f03_version(vmajor, vminor, vmicro)
  call xc_f90_info_name(xc_infox, xname)
  call xc_f90_info_name(xc_infoc, cname)
  write(6, '("# Libxc version: ",i1,".",i1,".",i1)') vmajor, vminor, vmicro
  write(6, '("# Exchange:    ", 1a)') trim(xname)
  write(6, '("# Correlation: ", 1a)') trim(cname)

end subroutine control_xc_init
!######################################################################
