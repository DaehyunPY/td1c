!######################################################################
subroutine hprod_trans(umat, wfn, twfn)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nact, nfun
  use mod_rad, only : nrad
  use mod_sph, only : lmax1

  implicit none
  complex(c_double_complex), intent(in) :: umat(1:nact, 1:nact)
  complex(c_double_complex), intent(in) :: wfn(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(out) :: twfn(1:(nrad-1), 0:lmax1, 1:nfun)
  integer(c_int) :: rdim, llr, ulr

  rdim = nrad - 1
  !$omp parallel default(shared) private(llr, ulr)
  !###########################
  call util_omp_disp(1, rdim, llr, ulr)
  call hprod_transp(umat, wfn, twfn, llr, ulr)
  !###########################
  !$omp end parallel

end subroutine hprod_trans
!######################################################################
subroutine hprod_transp(umat, wfn, twfn, llr, ulr)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : ncore, nact, nfun
  use mod_bas, only : mval
  use mod_rad, only : nrad
  use mod_sph, only : lmax1
  use mod_const, only : czero, ctwo

  implicit none
  integer(c_int), intent(in) :: llr, ulr
  complex(c_double_complex), intent(in) :: umat(1:nact, 1:nact)
  complex(c_double_complex), intent(in) :: wfn(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(out) :: twfn(1:(nrad-1), 0:lmax1, 1:nfun)

  integer(c_int) :: ifun, jfun, iact, jact, l, m, irad

  twfn(llr:ulr, 0:lmax1, 1:nfun) = czero

  ! umat-transformed orbitals
  do ifun = 1, ncore
     m = mval(ifun)
     do l = abs(m), lmax1
        do irad = llr, ulr
           twfn(irad, l, ifun) = twfn(irad, l, ifun) + wfn(irad, l, ifun) * ctwo
        end do
     end do
  end do

  do iact = 1, nact
     ifun = ncore + iact
     do jact = 1, nact
        jfun = ncore + jact
        if (mval(jfun) == mval(ifun)) then
           m = mval(ifun)
           do l = abs(m), lmax1
              do irad = llr, ulr
                 twfn(irad, l, ifun) = twfn(irad, l, ifun) &
                                      + wfn(irad, l, jfun) * umat(jact, iact)
              end do
           end do
        end if
     end do
  end do

end subroutine hprod_transp
!######################################################################
