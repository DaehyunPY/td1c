!///////////////////////////////////////////////////////////////////////
subroutine wfn_mask(orb)

  use, intrinsic :: iso_c_binding

  use mod_bas, only : mval
  use mod_ormas, only : nfcore, nfun
  use mod_sph, only : lmax1
  use mod_rad, only : nrad, xrad, rmask, mask

  implicit none
  complex(c_double_complex), intent(inout) :: orb(1:(nrad-1), 0:lmax1, 1:nfun)
  integer(c_int) :: ifun, irad, l, min_irad, llr, ulr

  ! determine min_irad
  min_irad = 1
  do irad = 1, nrad - 1
     if (xrad(irad) > rmask) then
        min_irad = irad
        exit
     end if
  end do

  !$omp parallel default(shared) private(llr, ulr)
  call util_omp_disp(min_irad, nrad - 1, llr, ulr)
  do ifun = nfcore + 1, nfun
     do l = abs(mval(ifun)), lmax1
        do irad = llr, ulr
           orb(irad, l, ifun) = orb(irad, l, ifun) * mask(irad)
        end do
     end do
  end do
  !$omp end parallel

end subroutine wfn_mask
!///////////////////////////////////////////////////////////////////////
