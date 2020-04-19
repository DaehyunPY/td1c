!#######################################################################
subroutine hprod_enepole(lfield, wfn, cic)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : nbas
  use mod_rad, only : nrad
  use mod_sph, only : lmax1
  use mod_control, only : igauge
  use mod_ormas, only : neltot, nfun, nfcore
  use mod_const, only : zero, two, czero, runit
  use mod_hprod, only : dip_exp, vel_exp, acc_exp
  use mod_hprod, only : orb, h0orb, h1orb, gorb, torb, den1

  implicit none
  real(c_double), intent(in) :: lfield(1:3, 1:3)
  complex(c_double_complex), intent(in) :: wfn(*)
  complex(c_double_complex), intent(in) :: cic(*)
  complex(c_double_complex), external :: hprod_trace1
  complex(c_double_complex) :: zfac
!
!  den1 -> d^p_a   ! 1rdm created in hprod_htot
!   orb ->    |p>  ! orbital copied in hprod_htot
!  gorb -> f_p|p>  ! g-fock created in hprod_htot
! h0orb -> z  |p>  ! dipole
! h1orb -> v_z|p>  ! dipole velocity
! 
  ! g-fock vector
  call hprod_gfock(lfield, wfn, cic)

  ! property vector
  call zclear_omp(nbas*nfun, h0orb)
  call zclear_omp(nbas*nfun, h1orb)
  call hprod_zprod_all(runit, orb, h0orb)
  call hprod_pzprod_all(runit, orb, h1orb)

  ! froze-core projection
  gorb(1:(nrad-1), 0:lmax1, 1:nfcore) = czero
  call hprod_projfc(.false., orb, gorb)

  if (igauge == 1) then
     zfac = -lfield(3, 2)
!bad, why?  torb(1:(nrad-1), 0:lmax1, 1:nfun) = zfac * h0orb(1:(nrad-1), 0:lmax1, 1:nfun)
     torb(1:(nrad-1), 0:lmax1, 1:nfun) = czero
     call zaxpy_omp(nbas*nfun, zfac, h0orb, torb)
     call hprod_projfc(.false., orb, torb)     
!good, why? gorb(1:(nrad-1), 0:lmax1, 1:nfcore) = torb(1:(nrad-1), 0:lmax1, 1:nfcore)
     call zcopy_omp(nbas*nfun, torb, gorb)

!bad, why?  torb(1:(nrad-1), 0:lmax1, 1:nfun) = - lfield(3, 2) * h0orb(1:(nrad-1), 0:lmax1, 1:nfun)
     torb(1:(nrad-1), 0:lmax1, 1:nfun) = czero
     call zaxpy_omp(nbas*nfun, zfac, h0orb, torb)
     call hprod_projfc(.true., orb, torb)
!bad, why?  gorb(1:(nrad-1), 0:lmax1, nfcore+1:nfun) = gorb(1:(nrad-1), 0:lmax1, nfcore+1:nfun) + torb(1:(nrad-1), 0:lmax1, nfcore+1:nfun)
     call zaxpy_omp(nbas*nfun, runit, torb, gorb)
  end if

  dip_exp = dble(hprod_trace1(.true., den1, orb, h0orb))
  vel_exp = two * aimag(hprod_trace1(.true., den1, h0orb, gorb))
  acc_exp = two * aimag(hprod_trace1(.true., den1, h1orb, gorb)) - neltot(3) * lfield(3, 2)

end subroutine hprod_enepole
!#######################################################################
subroutine hprod_enepole_old(lfield, wfn, cic)

  use, intrinsic :: iso_c_binding
  use mod_const, only : zero, two, runit
  use mod_ormas, only : neltot, nfun, nfcore
  use mod_bas, only : nbas
  use mod_control, only : igauge
  use mod_hprod, only : orb, h0orb, h1orb, gorb, den1
  use mod_hprod, only : dip_exp, vel_exp, acc_exp

  implicit none
  real(c_double), intent(in) :: lfield(1:3, 1:3)
  complex(c_double_complex), intent(in) :: wfn(*)
  complex(c_double_complex), intent(in) :: cic(*)
  complex(c_double_complex), external :: hprod_trace1
!
!  den1 -> d^p_a   ! 1rdm created in hprod_htot
!   orb ->    |p>  ! orbital copied in hprod_htot
!  gorb -> f_p|p>  ! g-fock created in hprod_htot
! h0orb -> z  |p>  ! dipole
! h1orb -> v_z|p>  ! dipole velocity
! 
  ! g-fock vector
  call hprod_gfock(lfield, wfn, cic)

  ! froze-core projection
  call hprod_projfc(.false., orb, gorb)

  ! property vector
  call zclear_omp(nbas*nfun, h0orb)
  call zclear_omp(nbas*nfun, h1orb)
  call hprod_zprod_dyn(runit, orb, h0orb)
  call hprod_pzprod_dyn(runit, orb, h1orb)

  ! dipole
  dip_exp = dble(hprod_trace1(.true., den1, orb, h0orb))

  ! dipole velocity
  vel_exp = two * aimag(hprod_trace1(.true., den1, h0orb, gorb))

  ! dipole acceleration
  if (igauge == 0) then
     acc_exp = two * aimag(hprod_trace1(.true., den1, h1orb, gorb))
  else
     if (nfcore > 0) then
        call hprod_projfc(.true., orb, h0orb)
        call zaxpy_omp(nbas*nfun, -lfield(3, 2), h0orb, gorb)
     end if
     acc_exp = two * aimag(hprod_trace1(.true., den1, h1orb, gorb)) &
             - lfield(3, 2) * (neltot(3) - 2 * nfcore)
  end if

end subroutine hprod_enepole_old
!#######################################################################
