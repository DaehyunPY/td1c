!################################################################################
subroutine ormas_ipx(max_ipx, ovlp, cic, ipx)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nfun, nstr_alph, nstr_beta, cic_old, nact, tdcc, nelact, neltot

  implicit none
  !--------------------------------------------------------------------
  integer(c_int), intent(in) :: max_ipx
  complex(c_double_complex), intent(in) :: ovlp(1:nfun, 1:nfun)
  complex(c_double_complex), intent(in) :: cic(1:*)
  real(c_double), intent(out) :: ipx(0:*)
  !--------------------------------------------------------------------

  if (tdcc) then
     ipx(0:max_ipx) = 0d0
     return
  end if

!  if (neltot(3) > 1 .and. &
!      (nstr_alph.ne.nstr_beta .or. &
!       nelact(1).ne.nelact(2))) then
!     stop 'ormas_ipx nyi for spin polarized case'
!  end if

  if (nact == 0) then
     call ormas_ipx_core(max_ipx, ovlp, ipx)
  else if (.not. cic_old) then
     call ormas_ipx_ras(max_ipx, ovlp, cic, ipx)
!     call ormas_denipx_ras(max_ipx, ovlp, cic, ipx)
  else
     call ormas_ipx_old(max_ipx, ovlp, cic, ipx)
  end if

end subroutine ormas_ipx
!################################################################################
