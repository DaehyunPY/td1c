!#######################################################################
subroutine hprod_gfock_dft(lfield, wfn, cic)

  use, intrinsic :: iso_c_binding
  use mod_const, only : runit
  use mod_bas, only : nbas, nbas2, ngrid
  use mod_control, only : igauge, PSP
  use mod_ormas, only : nfcore, neltot
  use mod_hprod, only : orb, orbg, h0orb, h1orb, gorb, gorbg, v2orb, v2orbg
  use mod_hprod, only : rho1, rho2, v2sph, v2ang, den1

  implicit none
  real(c_double), intent(in) :: lfield(1:3, 1:3)
  complex(c_double_complex), intent(in) :: wfn(*)
  complex(c_double_complex), intent(in) :: cic(*)
  complex(c_double_complex) :: zfield

  real(c_double) :: edftj, edftx
  complex(c_double_complex), allocatable :: rhoj(:)
  complex(c_double_complex), allocatable :: v2js(:)
  complex(c_double_complex), allocatable :: v2jg(:)
  real(c_double), allocatable :: rhoxc(:)

  allocate(rhoj(1:nbas2))
  allocate(v2js(1:nbas2))
  allocate(v2jg(1:ngrid))
  allocate(rhoxc(1:ngrid))

  call hprod_htot_init                     ! scratch initialization
  call ormas_mkden1(cic, den1)             ! 1RDM
  call hprod_orbin(lfield, wfn, orb, orbg) ! current orbitals

  ! 1e operators
  zfield = lfield(3, 1)
  call hprod_tprod_all(orb, h0orb)
  if (PSP) call hprod_projpp(runit, orb, h0orb)

  if (igauge == 0) call hprod_zprod_all(zfield, orb, h1orb)
  if (igauge == 1) call hprod_pzprod_all(zfield, orb, h1orb)

  ! meanfield operator
  call hprod_mkrhoxc(rhoj, rhoxc)
  call hprod_poisson1(rhoj, v2js)
  call bas_sph2ang2one(0, v2js, v2jg)
  call hprod_mfprod_dftj(v2jg, orbg, gorbg, edftj)
  call hprod_mfprod_dftx(rhoxc, orbg, gorbg, edftx)
  call bas_ang2sph1_dyn(gorbg, gorb);
  if (nfcore > 0) then
     call bas_ang2sph1_fc(gorbg, gorb);
  end if

  ! gfock vector
  call hprod_gfock_sum()

  deallocate(rhoxc)
  deallocate(v2jg)
  deallocate(v2js)
  deallocate(rhoj)

!DEBUG
!  use mod_rad, only : nrad
!  use mod_sph, only : lmax1, lmax2
!  integer(c_long) :: isph, irad, ifun
!  write(6, "('gorb and v2orb:')")
!  do isph = 0, lmax1
!     do irad = 1, nrad - 1
!        write(6, "(2i5)", advance = 'no') isph, irad
!        do ifun = 1, nfun
!           write(6, "(4f20.10)", advance = 'no') gorb(irad, isph, ifun), v2orb(irad, isph, ifun)
!        end do
!        write(6, *)
!     end do
!  end do
!  write(6, "('rho2:')")
!  do isph = 0, lmax2
!     do irad = 1, nrad - 1
!        write(6, "(2i5)", advance = 'no') isph, irad
!        write(6, "(8f20.10)", advance = 'no') rho2(irad, isph, 1,1), rho2(irad, isph, 1,2), rho2(irad, isph, 2,1), rho2(irad, isph, 2,2)
!        write(6, *)
!     end do
!  end do
!  write(6, "('v2sph:')")
!  do isph = 0, lmax2
!     do irad = 1, nrad - 1
!        write(6, "(2i5)", advance = 'no') isph, irad
!        write(6, "(8f20.10)", advance = 'no') v2sph(irad, isph, 1,1), v2sph(irad, isph, 1,2), v2sph(irad, isph, 2,1), v2sph(irad, isph, 2,2)
!        write(6, *)
!     end do
!  end do
!  stop
!DEBUG

end subroutine hprod_gfock_dft
!#######################################################################
