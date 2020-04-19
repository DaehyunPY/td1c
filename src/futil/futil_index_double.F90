!################################################################################
subroutine futil_index_double(dsc, n, vec, index)

  use, intrinsic :: iso_c_binding

  implicit none
  logical, intent(in) :: dsc
  integer(c_long), intent(in) :: n
  real(c_double), intent(in) :: vec(1:n)
  integer(c_long), intent(out) :: index(1:n)

  integer(c_long) :: i
  real(c_double) :: tmp
  real(c_double), parameter :: one = 1.d+0
  real(c_double), allocatable :: vec2(:)

!  n = size(vec)
!  n = ubound(vec)
!debug
!write(6,"('futil_index_double: n = ', i10)") n
!stop
!debug
  allocate(vec2(1:n))
  vec2(1:n) = vec(1:n)

  if (dsc) then
     tmp = minval(vec2) - one
     do i = 1, n
        index(i) = maxloc(vec2, dim = 1)
        vec2(index(i)) = tmp
     end do
  else
     tmp = maxval(vec2) + one
     do i = 1, n
        index(i) = minloc(vec2, dim = 1)
        vec2(index(i)) = tmp
     end do
  end if

  deallocate(vec2)

end subroutine futil_index_double
!################################################################################
