!#######################################################################
subroutine hprod_semicanonical(wfn)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : nbas
  use mod_control, only : exact3j, PSP
  use mod_const, only : zero, czero, runit, ctrue
  use mod_ormas, only : neltot, ncore, nact, nfun, nelact, norb_sub
  use mod_hprod, only : orb, orbg, h0orb, h1orb, gorb, gorbg, rho2, v2sph, int1e, int2e

  implicit none
  complex(c_double_complex), intent(inout) :: wfn(1:nbas,1:nfun)

  complex(c_double_complex) :: fock(1:nact,1:nact)
  complex(c_double_complex) :: umat(1:nact,1:nact)
  real(c_double) :: lfield(1:9) = zero
  integer(c_int) :: iact,jact,kact,ifun,jfun,ibas

  if (nact==0 .or. nelact(3)<=1) return

  call hprod_htot_init
  call hprod_orbin(lfield, wfn, orb, orbg)
  call hprod_tprod_dyn(orb, h0orb)
  if (PSP) call hprod_projpp(runit, lfield, orb, h0orb)
  if (neltot(3) >= 2) then
     if (.not.exact3j) then
        call hprod_mkrho2_dyn
        call hprod_mkv2mf_dyn
        call hprod_mfprod2_ene
        call bas_ang2sph1_dyn(gorbg, gorb);
     else
        call hprod_mkrho2_x3j_dyn
        call hprod_mkv2mf_x3j_dyn
        call hprod_mfprodx3j_ene
     end if
  end if
  call hprod_mkint1_sph(ctrue, orb, h0orb, h1orb, gorb)
  if (.not.exact3j) then
     call hprod_mkint2_sph(rho2, v2sph)
  else
     call hprod_mkint2_x3j(rho2, v2sph)
  end if

  do iact = 1, nact
     do jact = 1, nact
        fock(jact,iact) = int1e(jact,iact)
        do kact = 1, norb_sub(1)
           fock(jact,iact) = fock(jact,iact) + int2e(jact,iact,kact,kact)*2d0 &
                                             - int2e(jact,kact,kact,iact)
        end do
     end do
  end do

  call futil_diag_comp(.false., nact, fock, umat)
  wfn(1:nbas,(ncore+1):nfun) = 0d0
  do iact = 1, nact
     ifun = ncore + iact
     do jact = 1, nact
        jfun = ncore + jact
        call util_zaxpy(nbas, umat(jact,iact), orb(1,0,jfun), 1, wfn(1,ifun), 1)
     end do
  end do

end subroutine hprod_semicanonical
!#######################################################################
