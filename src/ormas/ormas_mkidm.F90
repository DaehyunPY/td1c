!################################################################################
subroutine ormas_mkidm(ovlp, cic, dens)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : cic_old, nfun, nstr_alph, nstr_beta

  implicit none
  !--------------------------------------------------------------------
  complex(c_double_complex), intent(in) :: ovlp(1:nfun, 1:nfun)
  complex(c_double_complex), intent(in) :: cic(1:nstr_alph, 1:nstr_beta)
  complex(c_double_complex), intent(out) :: dens(1:nfun, 1:nfun)
  !--------------------------------------------------------------------

  if (.not. cic_old) then
     stop 'ormas_mkidm_ras nyi.'
     call ormas_mkidm_ras(ovlp, cic, dens)
  else
     call ormas_mkidm_old(ovlp, cic, dens)
  end if

end subroutine ormas_mkidm
!################################################################################
