!######################################################################
subroutine futil_hforb_read(nval, lval, nrad, chirad)
!
  use, intrinsic :: iso_c_binding

! construct atomic hamiltonian in the fedvr basis
!
  implicit none
  integer(c_int), intent(in) :: nval, lval, nrad
  complex(c_double_complex), intent(inout) :: chirad(0:nrad)

  character(len = 9), parameter :: fname = "HFOrb.dat"
  integer(c_int), parameter :: ior = 20
  integer(c_int) :: nrad0, maxfun, lmax, l, ifun, irad
  complex(c_double_complex) :: trash
  complex(c_double_complex), parameter :: czero = (0.d+0, 0.d+0)
  complex(c_double_complex), parameter :: runit = (1.d+0, 0.d+0)
  complex(c_double_complex), parameter :: iunit = (0.d+0, 1.d+0)

  chirad(0:nrad) = czero

  open(unit = ior, file = trim(fname), status = 'old', form = 'formatted')
  rewind(unit = ior)

  read(ior, "(I25)") nrad0
  read(ior, "(I25)") maxfun
  read(ior, "(I25)") lmax
  if (nrad0 .ne. nrad) stop "nrad .ne. nrad0 @ futil_hforb_read."
  if (nval > maxfun) stop "n > maxfun @ futil_hforb_read."
  if (lval > lmax) stop "l > lmax @ futil_hforb_read."

  do l = 0, lval - 1
     do ifun = 1, maxfun
        do irad = 1, nrad - 1
           read(ior, "(2E25.15)") trash
        end do
     end do
  end do
  do ifun = 1, nval - 1
     do irad = 1, nrad - 1
        read(ior, "(2E25.15)") trash
     end do
  end do
  do irad = 1, nrad - 1
     read(ior, "(2E25.15)") chirad(irad)
  end do

  close(unit = ior)

end subroutine futil_hforb_read
!######################################################################
