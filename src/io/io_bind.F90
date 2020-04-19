!#######################################################################
subroutine io_bind(name_)

  use, intrinsic :: iso_c_binding
  use mod_control, only : name

  implicit none
  character(*), intent(in) :: name_

  name = name_

end subroutine io_bind
!#######################################################################
