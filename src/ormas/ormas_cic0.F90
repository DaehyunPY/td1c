!################################################################################
subroutine ormas_cic0(cic)

  use, intrinsic :: iso_c_binding
  use mod_control, only : name
  use mod_const, only : runit
  use mod_ormas, only : nact, lcic, ntot_alph_beta, cic_old
  use mod_ormas, only : map1to2_alph, map1to2_beta, tdcc

  implicit none
  complex(c_double_complex), intent(out) :: cic(1:*)

  if (nact == 0) return

  call zclear_omp(lcic, cic)
!!!  cic(1) = runit

  if (tdcc) then
     continue
  else if (cic_old) then
     cic(1) = runit
  else
!     cic(mapf_detx(map1to2_alph(1),map1to2_beta(1))) = runit
     cic(ntot_alph_beta(map1to2_beta(1))+map1to2_alph(1)) = runit
  end if

!DEBUG
!  write(6,"('alph = ',i10)") map1to2_alph(1)
!  write(6,"('beta = ',i10)") map1to2_beta(1)
!  write(6,"('idet = ',i10)") ntot_alph_beta(map1to2_beta(1))+map1to2_alph(1)
!  call ormas_cic_print(cic, trim(name)//".cicoeffg")
!  stop "STOP for debug @ ormas_cic0."
!DEBUG

end subroutine ormas_cic0
!################################################################################
