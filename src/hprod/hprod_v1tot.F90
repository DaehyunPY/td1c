!#######################################################################
subroutine hprod_v1tot(dtime, lfield, wfn, cic, hwfn, hcic)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : nbas
  use mod_hprod, only : int1e
  use mod_control, only : iprojfc, igauge, icomp
  use mod_const, only : zero, czero, runit, iunit
  use mod_ormas, only : nfcore, nfun, lcic, froz

  implicit none
  real(c_double), intent(in) :: dtime, lfield(1:3, 1:3)
  complex(c_double_complex), intent(in) :: wfn(1:nbas, 1:nfun)
  complex(c_double_complex), intent(in) :: cic(1:lcic)
  complex(c_double_complex), intent(inout) :: hwfn(1:nbas, 1:nfun)
  complex(c_double_complex), intent(inout) :: hcic(1:lcic)

  integer(c_int) :: ifun
  complex(c_double_complex) :: zfield, tfac

  zfield = lfield(3, 1)

  hwfn(1:nbas,1:nfun) = czero
  hcic(1:lcic) = czero

  if (igauge == 0) then
     call hprod_zprod_dyn(zfield, wfn, hwfn)
  else
     call hprod_pzprod_dyn(zfield, wfn, hwfn)
  end if

!  if (projhigh) then
!     call hprod_projhigh(hwfn)
!  end if

  ! sigma vector
  call hprod_v1tot_mkint1(wfn, hwfn, int1e)
  call ormas_hcic1(int1e, cic, hcic)

  if (icomp == 0) then
     tfac = - runit * dtime
  else
     tfac = - iunit * dtime
  end if
  do ifun = nfcore + 1, nfun
     if (froz(ifun) < 0) then
        hwfn(1:nbas,ifun) = czero
     else
        hwfn(1:nbas,ifun) = hwfn(1:nbas,ifun) * tfac
     end if
  end do
  hcic(1:lcic) = hcic(1:lcic) * tfac

end subroutine hprod_v1tot
!######################################################################
subroutine hprod_v1tot_mkint1(wfn, hwfn, int1e)

  use, intrinsic :: iso_c_binding
  use mod_const, only : czero
  use mod_bas, only : nbas, mval
  use mod_ormas, only : nact, ncore, nfun

  implicit none
  complex(c_double_complex), intent(in) :: wfn(1:nbas, 1:nfun)
  complex(c_double_complex), intent(in) :: hwfn(1:nbas, 1:nfun)
  complex(c_double_complex), intent(out) :: int1e(1:nact, 1:nact)

  integer(c_int), external :: util_omp_nproc
  integer(c_int), external :: util_omp_iproc
  complex(c_double_complex) :: tmp
  complex(c_double_complex), allocatable :: int1p(:,:,:)
  integer(c_int) :: nproc, iproc, iact, jact, ifun, jfun, ibas, llb, ulb

  nproc = util_omp_nproc()
  int1e(1:nact, 1:nact) = czero
  allocate(int1p(1:nact, 1:nact, 0:(nproc-1)))
  !$omp parallel default(shared) private(iproc, ifun, jfun, llb, ulb, tmp)
  !###########################
  iproc = util_omp_iproc()
  int1p(1:nact, 1:nact, iproc) = czero
  call util_omp_disp(1, nbas, llb, ulb)
  do iact = 1, nact
     ifun = ncore + iact
     do jact = 1, iact
        jfun = ncore + jact
        tmp = czero
        if (mval(jfun) == mval(ifun)) then
           ! <j|h0+h1+h2|i>
           do ibas = llb, ulb
              tmp = tmp + conjg(wfn(ibas, jfun)) * hwfn(ibas, ifun)
           end do
        end if
        int1p(jact, iact, iproc) = tmp
     end do
  end do
  !###########################
  !$omp end parallel
  do iproc = 0, nproc - 1
     do iact = 1, nact
        do jact = 1, iact
           int1e(jact, iact) = int1e(jact, iact) + int1p(jact, iact, iproc)
        end do
     end do
  end do
  do iact = 1, nact
     do jact = iact + 1, nact
        int1e(jact, iact) = conjg(int1e(iact, jact))
     end do
  end do
  deallocate(int1p)

end subroutine hprod_v1tot_mkint1
!#######################################################################
