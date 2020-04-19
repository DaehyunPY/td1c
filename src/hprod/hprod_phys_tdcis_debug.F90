!######################################################################
subroutine hprod_phys_tdcis_debug(dtime, radipx, lfield, wfn, cic, dipa, vela, acca)
  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad, xrad
  use mod_sph, only : lmax1
  use mod_bas, only : nbas
  use mod_ormas, only : nfun, nfcore, lcic
  use mod_hprod, only : dcic, h0orb, h1orb, orb0, orb
  use mod_hprod, only : ridm_tdcis, orb0rot, h1orb0
  use mod_hprod, only : h0orb0
  use mod_hprod, only : den1
  use mod_const, only : runit, iunit
  use mod_control, only : istdcis

  use mod_ormas, only : nfcore_tdcis, froz

  implicit none
  real(c_double), intent(in) :: dtime, radipx, lfield(1:3, 1:3)
  complex(c_double_complex), intent(in) :: wfn(1:*)
  complex(c_double_complex), intent(in) :: cic(1:lcic) !ci0
  real(c_double), intent(inout) :: dipa(1:3), vela(1:3), acca(1:3)

  complex(c_double_complex) :: ovlpin(1:nfun, 1:nfun)
  complex(c_double_complex) :: hwfn(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex) :: hcic(1:lcic)
  complex(c_double_complex) :: dotorb(1:(nrad-1), 0:lmax1, 1:nfun) 
  complex(c_double_complex) :: dotcic

  ! expected value
  complex(c_double_complex) :: dip, dip2, vel, vel2
  complex(c_double_complex) :: ip0, ip1
  integer(c_int) :: ifun
  integer(c_int) :: l, llr, ulr, irad

  ! write(6 , "('hprod_phys_tdcis_debug:')") 
  ! write(6 , *)dtime
  ! write(6 , *)radipx
  ! write(6 , *)lfield(3,2)
  ! write(6 , *)cic(1)
  ! write(6 , "('hprod_phys_tdcis_dip: debug1')") 
  call hprod_htot_tdcis(dtime, lfield, wfn, cic, hwfn, hcic)
  call hprod_mkovlp(xrad(nrad - 1), orb, orb, ridm_tdcis)
  write(6 , "('hprod_phys_tdcis_dip: debug3 ')") 

  !###########################
  ! call util_omp_disp(1, nrad - 1, llr, ulr)
  do ifun = nfcore_tdcis + 1, nfun
     ! if(froz(ifun)<0) cycle
     do l = 0, lmax1
        do irad = 1, nrad - 1
        ! do irad = llr, ulr
           dotorb(irad, l, ifun) = hwfn(irad, l, ifun)/dtime - iunit*(h0orb(irad, l, ifun) + h1orb(irad, l, ifun))
        end do
     end do
  end do
  !###########################

  stop

  call hprod_dipole_tdcis(lfield, dotorb, cic(1), dotcic, dip, dip2, vel, vel2)

  dotcic = dcic(1)/dtime

  dipa = 0d0
  vela = 0d0
  acca = 0d0
  ! dipa(1) = real(dip)
  ! vela(1) = real(dip2)
  ! vela(2) = real(vel)
  ! acca(1) = real(vel2)
  ! acca(2) = real(ip0)
  ! acca(3) = real(ip1)

end subroutine hprod_phys_tdcis_debug
!######################################################################
subroutine hprod_dotorb_tdcis_debug(dtime, hwfn, dotorb)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nfcore, nfun, froz
  use mod_rad, only : nrad
  use mod_sph, only : lmax1
  use mod_const, only : runit, iunit, zero


  use mod_hprod, only : dcic, h0orb, h1orb, orb0, orb
  use mod_ormas, only : nfcore_tdcis

  implicit none
  real(c_double), intent(in) :: dtime
  complex(c_double_complex), intent(in) :: hwfn(1:(nrad-1), 0:lmax1, 1:nfun)
  ! complex(c_double_complex), intent(in) :: h0orb(1:(nrad-1), 0:lmax1, 1:nfun)
  ! complex(c_double_complex), intent(in) :: h1orb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: dotorb(1:(nrad-1), 0:lmax1, 1:nfun)


  complex(c_double_complex) :: tmp
  integer(c_int) :: ifun, l, llr, ulr, irad

   ! dotorb = hwfn/dtime
  !$omp parallel default(shared) private(llr, ulr)
  !###########################
  ! call util_omp_disp(1, nrad - 1, llr, ulr)
  ! do ifun = nfcore_tdcis + 1, nfun
  !    if(froz(ifun)<0) cycle
  !    do l = 0, lmax1
  !       do irad = llr, ulr
  !          tmp = h0orb(irad, l, ifun) + h1orb(irad, l, ifun)
  !          dotorb(irad, l, ifun) = dotorb(irad, l, ifun) - iunit*tmp
  !       end do
  !    end do
  ! end do
  !###########################
  !$omp end parallel
end subroutine hprod_dotorb_tdcis_debug
!######################################################################
