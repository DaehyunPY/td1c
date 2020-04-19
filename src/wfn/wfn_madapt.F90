!///////////////////////////////////////////////////////////////////////
subroutine wfn_madapt(cic)

  use, intrinsic :: iso_c_binding

  use mod_bas, only : mtot
  use mod_const, only : czero
  use mod_ormas, only : mval_alph, mval_beta, nact
  use mod_ormas, only : lcic, ndetx, mapr_detx, tdcc

  implicit none
  complex(c_double_complex), intent(inout) :: cic(1:lcic)
  integer(c_int) :: istr, jstr, idet

  if (nact == 0) return
  !#### no longer needed
  !write(6, "('skip wfn_madapt.')")
  return
  !#### no longer needed

  do idet = 1, ndetx
     if (mval_alph(mapr_detx(1,idet)) + &
         mval_beta(mapr_detx(2,idet)) /= mtot) then
        cic(idet) = czero
     end if
  end do

  if (tdcc) then
     do idet = 1, ndetx
        if (mval_alph(mapr_detx(1,idet)) + &
            mval_beta(mapr_detx(2,idet)) /= mtot) then
           cic(ndetx+idet) = czero
        end if
     end do
  end if

end subroutine wfn_madapt
!///////////////////////////////////////////////////////////////////////
