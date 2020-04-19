!#######################################################################
subroutine hprod_dipole(lfield, wfn, cic, dip, vel, acc)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : nbas
  use mod_sph, only : lmax1
  use mod_control, only : 
  use mod_rad, only : nrad, nradfc, wrad
  use mod_ormas, only : neltot, nfun, nfcore
  use mod_const, only : zero, two, four, czero, runit
  use mod_control, only : igauge, fedvr_normalized
  use mod_hprod, only : den1, orb, orbg, h0orb, h1orb, gorb, torb, v2orb

  implicit none
  real(c_double), intent(in) :: lfield(1:3, 1:3)
  complex(c_double_complex), intent(in) :: wfn(*)
  complex(c_double_complex), intent(in) :: cic(*)
  complex(c_double_complex), external :: hprod_trace1
  real(c_double), intent(out) :: dip(1:3), vel(1:3), acc(1:3)

  integer(c_long) :: ifun, l, irad, llr, ulr
  real(c_double) :: vcorr1, vcorr2, acorr1, acorr2
  complex(c_double_complex) :: tmp, vcorr1p, acorr1p
!
!  den1 -> d^p_a   ! 1rdm created in hprod_htot
!   orb ->    |p>  ! orbital copied in hprod_gfock
!  gorb -> f_p|p>  ! g-fock created in hprod_gfock
! h0orb -> z  |p>  ! dipole
! h1orb -> v_z|p>  ! dipole velocity
!  torb            ! scratch
! v2orb            ! scratch
! 
  if (nfcore > 0) then
     ! g-fock vector
     if (.true.) then
        call hprod_gfock(lfield, wfn, cic)
     else
        ! ##### For DEBUG: no FC<=>Dyn eXchange #####
        call hprod_gfock_nofcx(lfield, wfn, cic)
        ! ###########################################
     end if
  else
     call ormas_mkden1(cic, den1);
     call hprod_orbin(lfield, wfn, orb, orbg)
  end if

  ! property vectors
  call zclear_omp(nbas*nfun, h0orb)
  call zclear_omp(nbas*nfun, h1orb)
  call zclear_omp(nbas*nfun, v2orb)
  call hprod_zprod_all(runit, orb, h0orb)
  call hprod_pzprod_all(runit, orb, h1orb)
  call hprod_azprod_all(runit, orb, v2orb)

  ! add iX term for velocity gauge case
  if (nfcore > 0 .and. igauge == 1) call zaxpy_omp(nbas*nfun, lfield(3, 2), h0orb, gorb)

  ! Ehrenfest formulus
  dip(1) = dble(hprod_trace1(.true., den1, orb, h0orb))
  vel(1) = dble(hprod_trace1(.true., den1, orb, h1orb))
  acc(1) = dble(hprod_trace1(.true., den1, orb, v2orb))
  if (igauge == 1) vel(1) = vel(1) + lfield(3, 3) * neltot(3)
  acc(1) = acc(1) - lfield(3, 2) * neltot(3)

  ! first correction term
  vcorr1 = zero
  acorr1 = zero
  if (nfcore > 0) then
     !$omp parallel default(shared) private(llr, ulr, tmp, vcorr1p, acorr1p) reduction(+:vcorr1, acorr1)
     !###########################
     call util_omp_disp(1, nradfc, llr, ulr)
     vcorr1p = czero
     acorr1p = czero
     if (fedvr_normalized) then
        do ifun = 1, nfcore
           do l = 0, lmax1
              do irad = llr, ulr
                 tmp = gorb(irad, l, ifun)
                 vcorr1p = vcorr1p + conjg(h0orb(irad, l, ifun)) * tmp
                 acorr1p = acorr1p + conjg(h1orb(irad, l, ifun)) * tmp
              end do
           end do
        end do
     else
        do ifun = 1, nfcore
           do l = 0, lmax1
              do irad = llr, ulr
                 tmp = gorb(irad, l, ifun) / wrad(irad)
                 vcorr1p = vcorr1p + conjg(h0orb(irad, l, ifun)) * tmp
                 acorr1p = acorr1p + conjg(h1orb(irad, l, ifun)) * tmp
              end do
           end do
        end do
     end if
     vcorr1 = vcorr1 - four * aimag(vcorr1p)
     acorr1 = acorr1 - four * aimag(acorr1p)
     !###########################
     !$omp end parallel
  end if

  ! second correction term
  vcorr2 = zero
  acorr2 = zero
  if (nfcore > 0) then
     !NEW
     !NEW  if (igauge == 1) call zaxpy_omp(nbas*nfun, lfield(3, 2), h0orb, gorb)
     !NEW
     call hprod_projfc(.true., orb, gorb)
     vcorr2 = - two * aimag(hprod_trace1(.true., den1, h0orb, gorb))
     acorr2 = - two * aimag(hprod_trace1(.true., den1, h1orb, gorb))
  end if

  dip(2) = dip(1)
  vel(2) = vel(1) + vcorr1
  acc(2) = acc(1) + acorr1

  dip(3) = dip(1)
  vel(3) = vel(1) + vcorr1 + vcorr2
  acc(3) = acc(1) + acorr1 + acorr2
!OLD
!OLD if (igauge == 1) acc_exp = acc_exp - two * nfcore * lfield(3, 2)
!OLD
!NEW
!NEW  if (igauge == 1) acc_exp = acc_exp + two * nfcore * lfield(3, 2)
!NEW

!DEBUG
!  write(6, "('hprod_dipole: ', 8F20.10)") vcorr1, vcorr2, vel(3), acorr1, acorr2, acc(3), acc(1), acc(2)
!DEBUG

end subroutine hprod_dipole
!#######################################################################
