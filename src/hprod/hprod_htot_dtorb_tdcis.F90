!######################################################################
! test FC use nfcore_tdcis instead of nfcore
!######################################################################
subroutine hprod_htot_dtorb_tdcis(dtime, orb, h0orb, h1orb, gorb, v2orb, cic, hwfn)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nfcore, nfun, froz
  use mod_bas, only : mval
  use mod_sph, only : lmax1
  use mod_rad, only : nrad, wrad, xrad, nradfc
  use mod_const, only : runit, iunit, zero
  use mod_control, only : fedvr_normalized, icomp, isplit
  use mod_hprod, only : xmat
  !tdcis
  use mod_ormas, only : lcic
  use mod_hprod, only : orb0rot, h1orb0, tdcis_eig
  use mod_control, only : tdcis_rvg
  use mod_ormas, only :nfcore_tdcis
  use mod_rad, only : nradgs
  !tdcis
  implicit none
  real(c_double), intent(in) :: dtime
  complex(c_double_complex), intent(in) ::   orb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: h0orb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: h1orb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) ::  gorb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: v2orb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: cic(1:lcic)
  complex(c_double_complex), intent(out) :: hwfn(1:(nrad-1), 0:lmax1, 1:nfun)

  complex(c_double_complex) :: tfac, h012, pmat, cmat
  integer(c_int) :: ifun, jfun, irad, l, llr, ulr
  complex(c_double_complex), allocatable :: invs(:,:)
  complex(c_double_complex), allocatable :: twfn(:,:,:)

  !tdcis
  complex(c_double_complex) :: gphichi(1:nfun,1:nfun)
  complex(c_double_complex) :: h0phichi(1:nfun,1:nfun)
  complex(c_double_complex) :: h1phiphi(1:nfun,1:nfun)
  complex(c_double_complex) :: h1phichi(1:nfun,1:nfun)
  !tdcis

  if (icomp == 0) then
     tfac = - runit * dtime
  else
     tfac = - iunit * dtime
  end if

  !projector not assuming orthonormality
  !write(6,"('hprod_htot_dtorb: P=|i>S^{-1}_{ij}<j|')")
  !allocate(invs(1:nfun,1:nfun))
  !allocate(twfn(1:(nrad-1), 0:lmax1, 1:nfun))
  !call hprod_mkovlp(xrad(nrad), orb, orb, invs)
  !call util_matinv_reg(nfun, zero, invs, invs)
  !twfn = zero
  !do ifun = 1, nfun
  !   do jfun = 1, nfun
  !      twfn(:,:,ifun) = twfn(:,:,ifun) + orb(:,:,jfun)*invs(jfun,ifun)
  !   end do
  !end do
  !projector not assuming orthonormality

  call hprod_mkovlp(xrad(nradgs), orb0rot, gorb, gphichi) ! <phi|G|chi>
  call hprod_mkovlp(xrad(nradgs), orb0rot, h0orb, h0phichi) ! <phi|h0|chi>
  call hprod_mkovlp(xrad(nradgs), orb0rot, h1orb0, h1phiphi) ! <phi|h1|phi>
  call hprod_mkovlp(xrad(nradgs), orb0rot, h1orb, h1phichi) ! <phi|h1|chi>

  !$omp parallel default(shared) private(h012, pmat, cmat, llr, ulr)
  !###########################
  call util_omp_disp(1, nrad - 1, llr, ulr)
  ! Q-contribution
  ! do ifun = nfcore + 1, nfun
  do ifun = nfcore_tdcis + 1, nfun
     if (froz(ifun) < 0) cycle
     if (isplit == 0) then
        ! <j|h0+h1+h2|i>
        do l = 0, lmax1
           do irad = llr, ulr
              h012 = h0orb(irad, l, ifun) + h1orb(irad, l, ifun) + gorb(irad, l, ifun) &
                   + h1orb0(irad, l, ifun)*sqrt(2d0)*cic(1)
              hwfn(irad, l, ifun) = h012 * tfac
           end do
        end do
     else if (isplit == 2) then
        write(6, '("isplit = 2 not yet implemented in tdcis")')
        stop 
        ! <j|h1+h2|i>
        do l = 0, lmax1
           do irad = llr, ulr
              h012 = v2orb(irad, l, ifun) + gorb(irad, l, ifun) + h1orb(irad, l, ifun)
              hwfn(irad, l, ifun) = h012 * tfac
           end do
        end do
     else if (isplit == 1) then
        ! <j|h2|i>
        do l = 0, lmax1
           do irad = llr, ulr
              h012 = gorb(irad, l, ifun) + h1orb0(irad, l, ifun)*sqrt(2d0)*cic(1)
              hwfn(irad, l, ifun) = h012 * tfac
           end do
        end do
     end if
  end do

  ! P-contribution
  ! do ifun = nfcore + 1, nfun
  do ifun = nfcore_tdcis + 1, nfun
     if (froz(ifun) < 0) cycle
     do jfun = 1 , nfun
        if (mval(ifun) == mval(jfun)) then
           pmat = h0phichi(jfun, ifun) + h1phichi(jfun, ifun) + gphichi(jfun, ifun) &
                + h1phiphi(jfun,ifun)*sqrt(2d0)*cic(1) 
           cmat = h1phiphi(jfun, ifun)
           do l = 0, lmax1
              do irad = llr, ulr
                 hwfn(irad, l, ifun) = hwfn(irad, l, ifun) - orb0rot(irad, l, jfun) * pmat * tfac - orb(irad, l, jfun) * cmat * tfac
                 !projector not assuming orthonormality
                 !hwfn(irad, l, ifun) = hwfn(irad, l, ifun) + twfn(irad, l, jfun) * xmat(jfun, ifun)
                 !projector not assuming orthonormality
              end do
           end do
        end if
     end do
  end do
  !###########################
  !$omp end parallel

  !projector not assuming orthonormality
  !deallocate(twfn)
  !deallocate(invs)
  !projector not assuming orthonormality

end subroutine hprod_htot_dtorb_tdcis
!######################################################################
subroutine hprod_htot_dtorb_tdcis_rvg(dtime, orb, h0orb, h1orb, gorb, v2orb, cic, afield, hwfn)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nfcore, nfun, froz
  use mod_bas, only : mval
  use mod_sph, only : lmax1
  use mod_rad, only : nrad, wrad, xrad, nradfc
  use mod_const, only : runit, iunit, zero
  use mod_control, only : fedvr_normalized, icomp, isplit
  use mod_hprod, only : xmat
  ! tdcis
  use mod_ormas, only : lcic
  use mod_hprod, only : orb0rot, h1orb0, tdcis_eig
  use mod_control, only : tdcis_rvg
  use mod_ormas, only : nfcore_tdcis
  use mod_rad, only : nradgs
  ! tdcis
  ! tdcis_rvg
  use mod_hprod, only : ezphi 
  ! tdcis_rvg

  implicit none
  real(c_double), intent(in) :: dtime
  complex(c_double_complex), intent(in) ::   orb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: h0orb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: h1orb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) ::  gorb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: v2orb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: cic(1:lcic)
  complex(c_double_complex), intent(in) :: afield
  complex(c_double_complex), intent(out) :: hwfn(1:(nrad-1), 0:lmax1, 1:nfun)

  complex(c_double_complex) :: tfac, h012, pmat, cmat
  integer(c_int) :: ifun, jfun, irad, l, llr, ulr
  complex(c_double_complex), allocatable :: invs(:,:)
  complex(c_double_complex), allocatable :: twfn(:,:,:)

  ! tdcis
  complex(c_double_complex) :: gphichi(1:nfun,1:nfun)
  complex(c_double_complex) :: h0phichi(1:nfun,1:nfun)
  complex(c_double_complex) :: h1phiphi(1:nfun,1:nfun)
  complex(c_double_complex) :: h1phichi(1:nfun,1:nfun)
  ! tdcis
  ! tdcis_rvg
  complex(c_double_complex) :: ezphichi(1:nfun,1:nfun)
  complex(c_double_complex) :: ezphiphi(1:nfun,1:nfun)
  ! tdcis_rvg

  if (icomp == 0) then
     tfac = - runit * dtime
  else
     tfac = - iunit * dtime
  end if

  !projector not assuming orthonormality
  !write(6,"('hprod_htot_dtorb: P=|i>S^{-1}_{ij}<j|')")
  !allocate(invs(1:nfun,1:nfun))
  !allocate(twfn(1:(nrad-1), 0:lmax1, 1:nfun))
  !call hprod_mkovlp(xrad(nrad), orb, orb, invs)
  !call util_matinv_reg(nfun, zero, invs, invs)
  !twfn = zero
  !do ifun = 1, nfun
  !   do jfun = 1, nfun
  !      twfn(:,:,ifun) = twfn(:,:,ifun) + orb(:,:,jfun)*invs(jfun,ifun)
  !   end do
  !end do
  !projector not assuming orthonormality

  ! tdcis_rvg
  call hprod_mkovlp(xrad(nradgs), ezphi, orb0rot,ezphiphi) ! E<phi|z|phi>
  call hprod_mkovlp(xrad(nradgs), ezphi, orb, ezphichi) ! E<phi|z|chi>
  ! tdcis_rvg
  call hprod_mkovlp(xrad(nradgs), orb0rot, gorb, gphichi) ! <phi|G|chi>
  call hprod_mkovlp(xrad(nradgs), orb0rot, h0orb, h0phichi) ! <phi|h0|chi>
  call hprod_mkovlp(xrad(nradgs), orb0rot, h1orb0, h1phiphi) ! <phi|h1|phi>
  call hprod_mkovlp(xrad(nradgs), orb0rot, h1orb, h1phichi) ! <phi|h1|chi>

  !$omp parallel default(shared) private(h012, pmat, cmat, llr, ulr)
  !###########################
  call util_omp_disp(1, nrad - 1, llr, ulr)
  ! Q-contribution
  ! do ifun = nfcore + 1, nfun
  do ifun = nfcore_tdcis + 1, nfun
     if (froz(ifun) < 0) cycle
     if (isplit == 0) then
        ! <j|h0+h1+h2|i>
        do l = 0, lmax1
           do irad = llr, ulr
              h012 = h0orb(irad, l, ifun) + h1orb(irad, l, ifun) + gorb(irad, l, ifun) &
                    + ezphi(irad, l, ifun)*sqrt(2d0)*cic(1) &
                    + abs(afield)*abs(afield)*orb(irad, l, ifun)*0.50
              hwfn(irad, l, ifun) = h012 * tfac
           end do
        end do
     else if (isplit == 2) then
        write(6, '("isplit = 2 not yet implemented in tdcis")')
        stop 
        ! <j|h1+h2|i>
        do l = 0, lmax1
           do irad = llr, ulr
              h012 = v2orb(irad, l, ifun) + gorb(irad, l, ifun) + h1orb(irad, l, ifun)
              hwfn(irad, l, ifun) = h012 * tfac
           end do
        end do
     else if (isplit == 1) then
        ! <j|h2|i>
        do l = 0, lmax1
           do irad = llr, ulr
              h012 = gorb(irad, l, ifun) + ezphi(irad, l, ifun)*sqrt(2d0)*cic(1) &
                   + abs(afield)*abs(afield)*orb(irad, l, ifun)*0.50
              hwfn(irad, l, ifun) = h012 * tfac
           end do
        end do
     end if
  end do

  ! P-contribution
  ! do ifun = nfcore + 1, nfun
  do ifun = nfcore_tdcis + 1, nfun
     if (froz(ifun) < 0) cycle
     do jfun = 1, nfun
        if (mval(ifun) == mval(jfun)) then
           pmat = h0phichi(jfun, ifun) + h1phichi(jfun, ifun) + gphichi(jfun, ifun) &
                + ezphiphi(jfun, ifun)*sqrt(2d0)*cic(1) + ezphichi(jfun, ifun)
           cmat = ezphiphi(jfun, ifun)
           do l = 0, lmax1
              do irad = llr, ulr
                 hwfn(irad, l, ifun) = hwfn(irad, l, ifun) - orb0rot(irad, l, jfun)*pmat*tfac - orb(irad, l, jfun)*cmat*tfac
                 !projector not assuming orthonormality
                 !hwfn(irad, l, ifun) = hwfn(irad, l, ifun) + twfn(irad, l, jfun) * xmat(jfun, ifun)
                 !projector not assuming orthonormality
              end do
           end do
        end if
     end do
  end do
  !###########################
  !$omp end parallel

  !projector not assuming orthonormality
  !deallocate(twfn)
  !deallocate(invs)
  !projector not assuming orthonormality

end subroutine hprod_htot_dtorb_tdcis_rvg
!######################################################################
