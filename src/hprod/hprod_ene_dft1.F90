!######################################################################
real(c_double) function hprod_ene_dft1(orb, horb)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : nbas
  use mod_ormas, only : nfcore, nocc
  use mod_const, only : two, czero

  implicit none
  complex(c_double_complex), intent(in) ::  orb(1:nbas, 1:*)
  complex(c_double_complex), intent(in) :: horb(1:nbas, 1:*)

  complex(c_double_complex) :: tmp, ene
  integer(c_long) :: ifun, ibas, llb, ulb

  ene = czero
  !$omp parallel default(shared) private(llb, ulb, tmp) reduction(+:ene)
  !###########################
  call util_omp_disp(1, nbas, llb, ulb)
  tmp = czero
  do ifun = nfcore + 1, nocc
     do ibas = llb, ulb
        tmp = tmp + conjg(orb(ibas, ifun)) * horb(ibas, ifun)
     end do
  end do
  ene = ene + tmp
  !###########################
  !$omp end parallel

  hprod_ene_dft1 = dble(ene) * two
  return

end function hprod_ene_dft1
!######################################################################
