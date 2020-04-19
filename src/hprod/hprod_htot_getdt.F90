!######################################################################
subroutine hprod_htot_getdt(dtime)
!
! final time derivatives are returned in 
! v2orb for orbitals, and
! dcic  for ci coefficients.
!
  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nfun
  use mod_const, only : runit
  use mod_bas, only : nbas, mval
  use mod_control, only : icomp, isplit
  use mod_hprod, only : orb, h0orb, h1orb, gorb, v2orb, dcic

  implicit none
  real(c_double), intent(in) :: dtime

  call hprod_htot_mkxmat_fc(dtime, orb, h0orb, h1orb, gorb, v2orb)
  call hprod_htot_mkxmat_cc(dtime, orb, h0orb, h1orb, gorb, v2orb)
  call hprod_htot_mkxmat_ca(dtime, orb, h0orb, h1orb, gorb, v2orb)
  call hprod_htot_mkxmat_aa(dtime, orb, h0orb, h1orb, gorb, v2orb)
  call hprod_htot_sum(dtime, orb, h0orb, h1orb, gorb, v2orb, dcic)

end subroutine hprod_htot_getdt
!######################################################################
subroutine hprod_htoto_getdt(dtime)
!
! final time derivatives are returned in 
! v2orb for orbitals.
!
  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nfun
  use mod_const, only : runit
  use mod_bas, only : nbas, mval
  use mod_control, only : icomp, isplit
  use mod_hprod, only : orb, h0orb, h1orb, gorb, v2orb

  implicit none
  real(c_double), intent(in) :: dtime

  call hprod_htot_mkxmat_fc(dtime, orb, h0orb, h1orb, gorb, v2orb)
  call hprod_htot_mkxmat_cc(dtime, orb, h0orb, h1orb, gorb, v2orb)
  call hprod_htot_mkxmat_ca(dtime, orb, h0orb, h1orb, gorb, v2orb)
  call hprod_htot_mkxmat_aa(dtime, orb, h0orb, h1orb, gorb, v2orb)
  call hprod_htoto_sum(dtime, orb, h0orb, h1orb, gorb, v2orb)

end subroutine hprod_htoto_getdt
!######################################################################
