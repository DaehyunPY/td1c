!#######################################################################
subroutine hprod_h1tot(dofc, lfield, wfn, hwfn)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : nbas
  use mod_const, only : runit
  use mod_ormas, only : nfun, nfcore, froz
  use mod_control, only : igauge, PSP

  implicit none
  logical(c_bool), intent(in) :: dofc
  real(c_double), intent(in) :: lfield(1:3, 1:3)
  complex(c_double_complex), intent(in) :: wfn(1:nbas,1:*)
  complex(c_double_complex), intent(inout) :: hwfn(1:nbas,1:*)
  
  logical(c_bool) :: dofcx
  integer(c_long) :: ifun
  complex(c_double_complex) :: zfield
  complex(c_double_complex), allocatable :: h1wfn(:,:)


  ! initialization
  allocate(h1wfn(1:nbas,1:nfun))
  zfield = lfield(3, 1)
  dofcx = dofc .and. nfcore > 0
  h1wfn = 0d0

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

  ! pseudopotential
  if (PSP) call hprod_projpp(runit, wfn, h1wfn)

  do ifun = nfcore + 1, nfun
     if (froz(ifun) < 0) cycle
     hwfn(1:nbas,ifun) = hwfn(1:nbas,ifun) + h1wfn(1:nbas,ifun)
  end do

end subroutine hprod_h1tot
!#######################################################################
