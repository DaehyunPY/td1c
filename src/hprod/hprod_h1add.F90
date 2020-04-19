!#######################################################################
subroutine hprod_h1add(dofc, zfac, lfield, wfn, hwfn)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : nbas
  use mod_const, only : czero
  use mod_control, only : igauge
  use mod_ormas, only : nfun, nfcore

  implicit none
  logical(c_bool), intent(in) :: dofc
  real(c_double), intent(in) :: lfield(1:3, 1:3)
  complex(c_double_complex), intent(in) :: zfac
  complex(c_double_complex), intent(in) :: wfn(1:*)
  complex(c_double_complex), intent(inout) :: hwfn(1:*)
  complex(c_double_complex), allocatable :: h1wfn(:)

  logical(c_bool) :: dofcx
  complex(c_double_complex) :: zfield

  allocate(h1wfn(1:nbas*nfun))
  h1wfn(1:nbas*nfun) = czero

  ! initialization
  zfield = lfield(3, 1)
  dofcx = dofc .and. nfcore > 0

  ! atomic hamiltonian
  call hprod_tprod_dyn(wfn, h1wfn)
  if (dofcx) call hprod_tprod_fc(wfn, h1wfn);

  ! laser hamiltonian
  if (igauge == 0) then
     call hprod_zprod_dyn(zfield, wfn, h1wfn)
     if (dofcx) call hprod_zprod_fc(zfield, wfn, h1wfn)
  else
     call hprod_pzprod_dyn(zfield, wfn, h1wfn)
     if (dofcx) call hprod_pzprod_fc(zfield, wfn, h1wfn)
  end if

  hwfn(1:nbas*nfun) = hwfn(1:nbas*nfun) + h1wfn(1:nbas*nfun) * zfac
  deallocate(h1wfn)

end subroutine hprod_h1add
!#######################################################################
