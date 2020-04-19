!#######################################################################
subroutine hprod_ipx(max_ipx, rad_ipx, lfield, wfn, cic, ipx)

  use, intrinsic :: iso_c_binding
  use mod_hprod, only : ovlp, orb, orbg

  implicit none
  integer(c_long), intent(in) :: max_ipx
  real(c_double), intent(in) :: rad_ipx
  real(c_double), intent(in) :: lfield(1:3, 1:3)
  complex(c_double_complex), intent(in) :: wfn(1:*)
  complex(c_double_complex), intent(in) :: cic(1:*)
  real(c_double), intent(out) :: ipx(1:*)

  call hprod_orbin(lfield, wfn, orb, orbg)
  call hprod_mkovlp(rad_ipx, orb, orb, ovlp);
  call ormas_ipx(max_ipx, ovlp, cic, ipx);

end subroutine hprod_ipx
!#######################################################################
subroutine hprod_ipd(max_ipx, rad_ipx, lfield, wfn, cic, ipx)

  use, intrinsic :: iso_c_binding
  use mod_hprod, only : ovlp, orb, orbg

  implicit none
  integer(c_long), intent(in) :: max_ipx
  real(c_double), intent(in) :: rad_ipx
  real(c_double), intent(in) :: lfield(1:3, 1:3)
  complex(c_double_complex), intent(in) :: wfn(1:*)
  complex(c_double_complex), intent(in) :: cic(1:*)
  real(c_double), intent(out) :: ipx(1:*)

  call hprod_orbin(lfield, wfn, orb, orbg)
  call hprod_mkovlp(rad_ipx, orb, orb, ovlp);
  call ormas_ipd(max_ipx, ovlp, cic, ipx);

end subroutine hprod_ipd
!#######################################################################
