!#######################################################################
subroutine hprod_v1ext(dofc, zfac, lfield, wfn, hwfn)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : nbas
  use mod_control, only : igauge
  use mod_ormas, only : nfun, nfcore

  implicit none
  logical(c_bool), intent(in) :: dofc
  complex(c_double_complex), intent(in) :: zfac
  real(c_double), intent(in) :: lfield(1:3, 1:3)
  complex(c_double_complex), intent(in) :: wfn(*)
  complex(c_double_complex), intent(inout) :: hwfn(*)

  logical(c_bool) :: dofcx
  complex(c_double_complex) :: zfield

  ! initialization
  zfield = lfield(3, 1) * zfac
  dofcx = dofc .and. nfcore > 0

  ! laser hamiltonian
  if (igauge == 0) then
     call hprod_zprod_dyn(zfield, wfn, hwfn)
     if (dofcx) call hprod_zprod_fc(zfield, wfn, hwfn)
  else
     call hprod_pzprod_dyn(zfield, wfn, hwfn)
     if (dofcx) call hprod_pzprod_fc(zfield, wfn, hwfn)
  end if

end subroutine hprod_v1ext
!######################################################################
