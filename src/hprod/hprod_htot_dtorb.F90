!######################################################################
subroutine hprod_htot_dtorb(dtime, orb, h0orb, h1orb, gorb, v2orb, hwfn)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nfcore, nfun, froz
  use mod_bas, only : mval
  use mod_sph, only : lmax1
  use mod_rad, only : nrad, wrad, xrad
  use mod_const, only : runit, iunit, zero
  use mod_control, only : fedvr_normalized, icomp, isplit
  use mod_hprod, only : xmat

  implicit none
  real(c_double), intent(in) :: dtime
  complex(c_double_complex), intent(in) ::   orb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: h0orb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: h1orb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) ::  gorb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: v2orb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(out) :: hwfn(1:(nrad-1), 0:lmax1, 1:nfun)

  complex(c_double_complex) :: tfac, h012
  integer(c_long) :: ifun, jfun, irad, l, llr, ulr
  complex(c_double_complex), allocatable :: invs(:,:)
  complex(c_double_complex), allocatable :: twfn(:,:,:)

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

  !$omp parallel default(shared) private(h012, llr, ulr)
  !###########################
  call util_omp_disp(1, nrad - 1, llr, ulr)
  ! Q-contribution
  do ifun = nfcore + 1, nfun
     if (froz(ifun) < 0) cycle
     if (isplit == 0) then
        ! <j|h0+h1+h2|i>
        do l = 0, lmax1
           do irad = llr, ulr
              h012 = v2orb(irad, l, ifun) + gorb(irad, l, ifun) + h1orb(irad, l, ifun) + h0orb(irad, l, ifun)
              hwfn(irad, l, ifun) = h012 * tfac
           end do
        end do
     else if (isplit == 2) then
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
              h012 = v2orb(irad, l, ifun) + gorb(irad, l, ifun)
              hwfn(irad, l, ifun) = h012 * tfac
           end do
        end do
     end if
  end do

  ! P-contribution
  do ifun = nfcore + 1, nfun
     if (froz(ifun) < 0) cycle
     do jfun = 1, nfun
        if (mval(ifun) == mval(jfun)) then
           do l = 0, lmax1
              do irad = llr, ulr
                 hwfn(irad, l, ifun) = hwfn(irad, l, ifun) + orb(irad, l, jfun) * xmat(jfun, ifun)
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

end subroutine hprod_htot_dtorb
!######################################################################
