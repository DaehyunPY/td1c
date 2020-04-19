!################################################################################
subroutine ormas_rm_singles(cic)

  use, intrinsic :: iso_c_binding

  use mod_ormas, only : nact, ntot_alph_beta, lcic
  use mod_ormas, only : nstr_beta, dist_str_alph, dist_str_beta
  use mod_ormas, only : llstr_alph_beta, nstr_alph_beta, ormas_donly

  implicit none
  complex(c_double_complex), intent(inout) :: cic(1:lcic)
  integer(c_long) :: istr, jstr, ifun, jfun

  if (nact == 0) return
  if (.not. ormas_donly) return

  !$omp parallel default(shared) private(istr,jstr)
  !$omp do
  do istr = 1, nstr_beta
     do jstr = llstr_alph_beta(istr), llstr_alph_beta(istr)+nstr_alph_beta(istr)-1
        if (dist_str_alph(1,jstr)+dist_str_beta(1,istr)-2 == 1) cic(ntot_alph_beta(istr)+jstr) = 0d0
!        write(6, "('ormas_rm_singles',3i7,2f20.10)") &
!             dist_str_alph(1,jstr)-1,dist_str_beta(1,istr)-1, &
!             dist_str_alph(1,jstr) + dist_str_beta(1,istr)-2, &
!             cic(ntot_alph_beta(istr)+jstr)
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ormas_rm_singles
!################################################################################
