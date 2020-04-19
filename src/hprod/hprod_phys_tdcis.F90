!######################################################################
subroutine hprod_phys_tdcis(dtime, radipx, lfield, wfn, cic, dipa, vela, acca)
  ! dipole, velocity, momenta, ipx
  ! called every time nstep_dip
  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad, xrad, ecs_flag, irad_ecs
  use mod_sph, only : lmax1
  use mod_bas, only : nbas
  use mod_ormas, only : nfun, nfcore, lcic
  use mod_hprod, only : dcic, h0orb, h1orb, orb0, orb
  use mod_hprod, only : ridm_tdcis, orb0rot, h1orb0
  use mod_hprod, only : h0orb0
  use mod_hprod, only : den1
  use mod_const, only : runit, iunit
  use mod_control, only : istdcis

  use mod_ormas, only : nfcore_tdcis

  implicit none
  real(c_double), intent(in) :: dtime, radipx, lfield(1:3, 1:3)
  complex(c_double_complex), intent(in) :: wfn(1:*)
  complex(c_double_complex), intent(in) :: cic(1:lcic) !ci0
  real(c_double), intent(out) :: dipa(1:3), vela(1:3), acca(1:3)

  ! expected value
  complex(c_double_complex) :: dip, dip2, vel, vel2
  complex(c_double_complex) :: ip0, ip1
  integer(c_int) :: ifun

  complex(c_double_complex) :: dotcic
!td1c_sato
!  complex(c_double_complex) :: ovlpin(1:nfun, 1:nfun)
!  complex(c_double_complex) :: hwfn(1:nbas*nfun)
!  complex(c_double_complex) :: hcic(1:lcic)
!  complex(c_double_complex) :: dotorb(1:(nrad - 1), 0:lmax1, 1:nfun) 
  real(c_double) :: radmax
  complex(c_double_complex), allocatable :: ovlpin(:,:)
  complex(c_double_complex), allocatable :: hwfn(:)
  complex(c_double_complex), allocatable :: hcic(:)
  complex(c_double_complex), allocatable :: dotorb(:,:,:) 
  allocate(ovlpin(1:nfun, 1:nfun))
  allocate(hwfn(1:nbas*nfun))
  allocate(hcic(1:lcic))
  allocate(dotorb(1:(nrad - 1), 0:lmax1, 1:nfun) )
!td1c_sato


  !!!debug
  ! call hprod_mkovlp(xrad(nrad), orb, orb, ridm_tdcis)
  ! write(6 , "('hprod_phys_tdcis_dip: debug ovlp')") 
  ! do ifun = 1, nfun
  !    write(6 , "(2es20.10)") ridm_tdcis(ifun, ifun)
  ! end do
  ! return
  !!!debug

!write(6,"('for debug F0.')")
!stop

  ! dotorb
  call hprod_htot_tdcis(dtime, lfield, wfn, cic, hwfn, hcic)

!write(6,"('for debug F1.')")
!stop


!tdcis_sato
  if (ecs_flag == 0) then
     radmax = xrad(nrad)
  else 
     radmax = xrad(irad_ecs-1)
  end if
!tdcis_sato

  dotorb = 0d0
  call hprod_mkovlp(radmax, orb, orb, ridm_tdcis)
!write(6,"('for debug F2.')")
!stop
  call hprod_dotorb_tdcis(dtime, hwfn, h0orb, h1orb, dotorb)
!write(6,"('for debug F3.')")
!stop
  dotcic = dcic(1)/dtime
  ! dotorb
 
  ! dipole
  call hprod_dipole_tdcis(lfield, dotorb, cic(1), dotcic, dip, dip2, vel, vel2)
  ! dipole
!write(6,"('for debug F4.')")
!stop
 
  !ipx
  ip0 = 0d0
  ip1 = 0d0
  call hprod_mkovlp(radipx, orb, orb, ovlpin)
  ! do ifun = nfcore + 1, nfun
  do ifun = nfcore_tdcis + 1, nfun
     ip0 = ip0 + ridm_tdcis(ifun, ifun)
     ip1 = ip1 + ovlpin(ifun, ifun)
  end do
  !ipx
 
  ! debug
  ! write(6 , "('hprod_phys_tdcis_dip: debug ')") 
  ! do ifun = 1, nfun
  !    write(6 , "(4es20.10)") ridm_tdcis(ifun, ifun), ovlpin(ifun, ifun)
  ! end do
  ! debug


  ! debug
  ! write(6 , "('hprod_phys_tdcis_ipx:', 6es20.10)") cic(1), ridm_tdcis(1,1), ip1
  ! write(6 , "('hprod_phys_tdcis_dip:', 8es20.10)") dip, dip2, vel, vel2
  ! debug

  dipa(1) = real(dip)
  vela(1) = real(dip2)
  vela(2) = real(vel)
  acca(1) = real(vel2)
  acca(2) = real(ip0)
  acca(3) = real(ip1)

!tdcis_sato
  deallocate(ovlpin)
  deallocate(hwfn)
  deallocate(hcic)
  deallocate(dotorb)
!tdcis_sato

!write(6,"('for debug F5.')")
!stop

end subroutine hprod_phys_tdcis
!######################################################################
subroutine hprod_dotorb_tdcis(dtime, hwfn, h0orb, h1orb, dotorb)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nfcore, nfun, froz
  use mod_rad, only : nrad
  use mod_sph, only : lmax1
  use mod_const, only : runit, iunit, zero

  use mod_ormas, only : nfcore_tdcis

  implicit none
  real(c_double), intent(in) :: dtime
  complex(c_double_complex), intent(in) :: hwfn(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: h0orb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: h1orb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(inout) :: dotorb(1:(nrad-1), 0:lmax1, 1:nfun)

  complex(c_double_complex) :: tmp
  integer(c_int) :: ifun, l, llr, ulr, irad

  dotorb = hwfn/dtime
  !$omp parallel default(shared) private(llr, ulr)
  !###########################
  call util_omp_disp(1, nrad - 1, llr, ulr)
  do ifun = nfcore_tdcis + 1, nfun
     if(froz(ifun)<0) cycle
     do l = 0, lmax1
        do irad = llr, ulr
           tmp = h0orb(irad, l, ifun) + h1orb(irad, l, ifun)
           dotorb(irad, l, ifun) = dotorb(irad, l, ifun) - iunit*tmp
        end do
     end do
  end do
  !###########################
  !$omp end parallel

end subroutine hprod_dotorb_tdcis
!######################################################################
