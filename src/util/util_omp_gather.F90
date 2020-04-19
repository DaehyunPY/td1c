!###########################################################
subroutine util_omp_gather(ndat, datp, dat)

  use, intrinsic :: iso_c_binding
  use mod_const, only : czero

  implicit none
  integer(c_long), intent(in) :: ndat
  complex(c_double_complex), intent(in) :: datp(1:ndat, 0:*)
  complex(c_double_complex), intent(out) :: dat(1:ndat)
  integer(c_long), external :: util_omp_nproc
  integer(c_long) :: idat, iproc, nproc

  nproc = util_omp_nproc()

  dat(1:ndat) = czero
  do iproc = 0, nproc - 1
     do idat = 1, ndat
        dat(idat) = dat(idat) + datp(idat, iproc)
     end do
  end do

end subroutine util_omp_gather
!###########################################################
