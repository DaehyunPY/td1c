!######################################################################
subroutine hprod_htot_sum(dtime, orb, h0orb, h1orb, gorb, v2orb, dcic)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : mval
  use mod_sph, only : lmax1
  use mod_rad, only : nrad, wrad
  use mod_const, only : runit, iunit
  use mod_control, only : fedvr_normalized, icomp, isplit
  use mod_hprod, only : xmat
  use mod_ormas, only : nfun, nfcore, lcic

  implicit none
  real(c_double), intent(in) :: dtime
  complex(c_double_complex), intent(in) ::   orb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: h0orb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: h1orb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(inout) ::  gorb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(inout) :: v2orb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(inout) :: dcic(1:lcic)

  complex(c_double_complex) :: tfac
  integer(c_long) :: idet, ifun, jfun, irad, l, llr, ulr

  if (icomp == 0) then
     tfac = - runit * dtime
  else
     tfac = - iunit * dtime
  end if

  !$omp parallel default(shared)
  !$omp do
  !###########################
  do idet = 1, lcic
     dcic(idet) = tfac * dcic(idet)
  end do
  !###########################
  !$omp end do
  !$omp end parallel

  !$omp parallel default(shared) private(llr, ulr)
  !###########################
  call util_omp_disp(1, nrad - 1, llr, ulr)
  ! Q-contribution
  do ifun = 1, nfun
     if (isplit == 0) then
        ! <j|h0+h1+h2|i>
        do l = 0, lmax1
           do irad = llr, ulr
              gorb(irad, l, ifun) = v2orb(irad, l, ifun) + gorb(irad, l, ifun) + h1orb(irad, l, ifun) + h0orb(irad, l, ifun)
           end do
        end do
     else if (isplit == 2) then
        ! <j|h1+h2|i>
        do l = 0, lmax1
           do irad = llr, ulr
              gorb(irad, l, ifun) = v2orb(irad, l, ifun) + gorb(irad, l, ifun) + h1orb(irad, l, ifun)
           end do
        end do
     else if (isplit == 1) then
        ! <j|h2|i>
        do l = 0, lmax1
           do irad = llr, ulr
              gorb(irad, l, ifun) = v2orb(irad, l, ifun) + gorb(irad, l, ifun)
           end do
        end do
     end if
  end do

  ! P-contribution
  do ifun = nfcore + 1, nfun
     if (fedvr_normalized) then
        do l = 0, lmax1
           do irad = llr, ulr
              v2orb(irad, l, ifun) = gorb(irad, l, ifun) * tfac
           end do
        end do
     else
        do l = 0, lmax1
           do irad = llr, ulr
              v2orb(irad, l, ifun) = gorb(irad, l, ifun) * tfac / wrad(irad)
           end do
        end do
     end if
     do jfun = 1, nfun
        if (mval(ifun) == mval(jfun)) then
           do l = 0, lmax1
              do irad = llr, ulr
                 v2orb(irad, l, ifun) = v2orb(irad, l, ifun) + orb(irad, l, jfun) * xmat(jfun, ifun)
              end do
           end do
        end if
     end do
  end do
  !###########################
  !$omp end parallel

end subroutine hprod_htot_sum
!######################################################################
subroutine hprod_htoto_sum(dtime, orb, h0orb, h1orb, gorb, v2orb)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : mval
  use mod_sph, only : lmax1
  use mod_rad, only : nrad, wrad
  use mod_const, only : runit, iunit
  use mod_control, only : fedvr_normalized, icomp, isplit
  use mod_hprod, only : xmat
  use mod_ormas, only : nfun, nfcore

  implicit none
  real(c_double), intent(in) :: dtime
  complex(c_double_complex), intent(in) ::   orb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: h0orb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: h1orb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(inout) ::  gorb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(inout) :: v2orb(1:(nrad-1), 0:lmax1, 1:nfun)

  complex(c_double_complex) :: tfac
  integer(c_long) :: ifun, jfun, irad, l, llr, ulr

  if (icomp == 0) then
     tfac = - runit * dtime
  else
     tfac = - iunit * dtime
  end if

  !$omp parallel default(shared) private(llr, ulr)
  !###########################
  call util_omp_disp(1, nrad - 1, llr, ulr)
  ! Q-contribution
  do ifun = 1, nfun
     if (isplit == 0) then
        ! <j|h0+h1+h2|i>
        do l = 0, lmax1
           do irad = llr, ulr
              gorb(irad, l, ifun) = v2orb(irad, l, ifun) + gorb(irad, l, ifun) + h1orb(irad, l, ifun) + h0orb(irad, l, ifun)
           end do
        end do
     else if (isplit == 2) then
        ! <j|h1+h2|i>
        do l = 0, lmax1
           do irad = llr, ulr
              gorb(irad, l, ifun) = v2orb(irad, l, ifun) + gorb(irad, l, ifun) + h1orb(irad, l, ifun)
           end do
        end do
     else if (isplit == 1) then
        ! <j|h2|i>
        do l = 0, lmax1
           do irad = llr, ulr
              gorb(irad, l, ifun) = v2orb(irad, l, ifun) + gorb(irad, l, ifun)
           end do
        end do
     end if
  end do

  ! P-contribution
  do ifun = nfcore + 1, nfun
     if (fedvr_normalized) then
        do l = 0, lmax1
           do irad = llr, ulr
              v2orb(irad, l, ifun) = gorb(irad, l, ifun) * tfac
           end do
        end do
     else
        do l = 0, lmax1
           do irad = llr, ulr
              v2orb(irad, l, ifun) = gorb(irad, l, ifun) * tfac / wrad(irad)
           end do
        end do
     end if
     do jfun = 1, nfun
        if (mval(ifun) == mval(jfun)) then
           do l = 0, lmax1
              do irad = llr, ulr
                 v2orb(irad, l, ifun) = v2orb(irad, l, ifun) + orb(irad, l, jfun) * xmat(jfun, ifun)
              end do
           end do
        end if
     end do
  end do
  !###########################
  !$omp end parallel

end subroutine hprod_htoto_sum
!######################################################################
subroutine hprod_htot_sum0(dtime, orb, h0orb, h1orb, gorb, v2orb, dcic)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : mval
  use mod_sph, only : lmax1
  use mod_rad, only : nrad, wrad
  use mod_const, only : runit, iunit
  use mod_control, only : fedvr_normalized, icomp
  use mod_hprod, only : xmat
  use mod_ormas, only : nfun, nfcore, lcic

  implicit none
  real(c_double), intent(in) :: dtime
  complex(c_double_complex), intent(in) ::   orb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: h0orb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: h1orb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(inout) ::  gorb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(inout) :: v2orb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(inout) :: dcic(1:lcic)

  complex(c_double_complex) :: tfac
  integer(c_long) :: idet, ifun, jfun, irad, l, llr, ulr

  if (icomp == 0) then
     tfac = - runit * dtime
  else
     tfac = - iunit * dtime
  end if

  !$omp parallel default(shared)
  !$omp do
  !###########################
  do idet = 1, lcic
     dcic(idet) = tfac * dcic(idet)
  end do
  !###########################
  !$omp end do
  !$omp end parallel

  !$omp parallel default(shared) private(llr, ulr)
  !###########################
  call util_omp_disp(1, nrad - 1, llr, ulr)
  ! Q-contribution
  do ifun = 1, nfun
     do l = 0, lmax1
        do irad = llr, ulr
           gorb(irad, l, ifun) = v2orb(irad, l, ifun) + gorb(irad, l, ifun) + h1orb(irad, l, ifun) + h0orb(irad, l, ifun)
        end do
     end do
  end do

  ! P-contribution
  do ifun = nfcore + 1, nfun
     if (fedvr_normalized) then
        do l = 0, lmax1
           do irad = llr, ulr
              v2orb(irad, l, ifun) = gorb(irad, l, ifun) * tfac
           end do
        end do
     else
        do l = 0, lmax1
           do irad = llr, ulr
              v2orb(irad, l, ifun) = gorb(irad, l, ifun) * tfac / wrad(irad)
           end do
        end do
     end if
     do jfun = 1, nfun
        if (mval(ifun) == mval(jfun)) then
           do l = 0, lmax1
              do irad = llr, ulr
                 v2orb(irad, l, ifun) = v2orb(irad, l, ifun) + orb(irad, l, jfun) * xmat(jfun, ifun)
              end do
           end do
        end if
     end do
  end do
  !###########################
  !$omp end parallel

end subroutine hprod_htot_sum0
!######################################################################
