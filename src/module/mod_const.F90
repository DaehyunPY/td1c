!################################################################################
module mod_const

  use, intrinsic :: iso_c_binding

  implicit none

  real(c_double), parameter :: zero  = 0.d+0
  real(c_double), parameter :: one   = 1.d+0
  real(c_double), parameter :: two   = 2.d+0
  real(c_double), parameter :: three = 3.d+0
  real(c_double), parameter :: four  = 4.d+0
  real(c_double), parameter :: five  = 5.d+0
  real(c_double), parameter :: six   = 6.d+0
  real(c_double), parameter :: eight = 8.d+0
  real(c_double), parameter :: half  = 5.d-1
  real(c_double), parameter :: third = 1.d+0 / 3.d+0
  real(c_double), parameter :: quart = 1.d+0 / 4.d+0
  real(c_double), parameter :: tenth = 1.d-1
  real(c_double), parameter :: pi = 3.14159265358979323846264338327950288d+0
  real(c_double), parameter :: au2fs = 0.02418884254d+0
  real(c_double), parameter :: fac_lda = (three / pi) ** third
  complex(c_double_complex), parameter :: czero = (zero, zero)
  complex(c_double_complex), parameter :: runit = (one, zero)
  complex(c_double_complex), parameter :: iunit = (zero, one)
  complex(c_double_complex), parameter :: chalf = (half, zero)
  complex(c_double_complex), parameter :: ctwo  = (two, zero)
  complex(c_double_complex), parameter :: cfour = (four, zero)

  logical(c_bool), parameter :: ctrue = .true.
  logical(c_bool), parameter :: cfalse = .false.
  real(c_double), parameter :: thrwfn = 1.D-15
  real(c_double), parameter :: thrdav = 1.D-12
  integer(c_long), parameter :: cmf_maxcic = 20
  integer(c_long), parameter :: maxproc = 256

end module mod_const
!################################################################################
