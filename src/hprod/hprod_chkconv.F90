!#######################################################################
subroutine hprod_chkconv(lfield, wfn, cic, ene, grad)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : nbas
  use mod_rad, only : nrad,xrad
  use mod_const, only : zero,one,two,cfalse
  use mod_hprod, only : ene_tot,den1,orb,gorb,ovlp
  use mod_ormas, only : nfcore,ncore,nact,nfun,nocc,lcic

  implicit none
  real(c_double), intent(in) :: lfield(1:3, 1:3)
  complex(c_double_complex), intent(in) :: wfn(1:*)
  complex(c_double_complex), intent(in) :: cic(1:*)
  real(c_double), intent(out) :: ene, grad(1:*)
  complex(c_double_complex), allocatable :: hwfn(:)
  complex(c_double_complex), allocatable :: hcic(:)
  complex(c_double_complex) :: cval
  integer(c_long) :: ifun,iact,jact

  allocate(hwfn(1:nbas*nfun))
  allocate(hcic(1:lcic))

!debug  call hprod_energy(lfield, wfn, cic)
!debug  ene = ene_tot

  call hprod_htot(one, lfield, wfn, cic, hwfn, hcic)
  ene = ene_tot
  call zdotc_omp(lcic, hcic, hcic, cval)
  grad(1) = sqrt(dble(cval)) / lcic

  call hprod_gfock(lfield, wfn, cic)
  call hprod_projall(cfalse, orb, gorb)
  call hprod_mkovlp(xrad(nrad), gorb, gorb, ovlp)

  cval = zero
  do ifun = nfcore + 1, nocc
     cval = cval + ovlp(ifun,ifun)
  end do
  grad(2) = sqrt(dble(cval)) / nfun

  cval = zero
  do ifun = nfcore + 1, ncore
     cval = cval + two * ovlp(ifun,ifun)
  end do
  do iact = 1, nact
     do jact = 1, nact
        cval = cval + ovlp(ncore+iact,ncore+jact) * den1(jact,iact)
     end do
  end do
  grad(3) = sqrt(dble(cval))

  !DEBUG
!  do ifun = 1, nfun
!     write(6, "('diag_gfock:', i5, 2e25.15)") ifun, ovlp(ifun,ifun)
!  end do
  !DEBUG

  deallocate(hwfn)
  deallocate(hcic)

end subroutine hprod_chkconv
!#######################################################################
