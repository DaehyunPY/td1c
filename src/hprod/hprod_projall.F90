!######################################################################
subroutine hprod_projall(inner, wfn, hwfn)
!
! for normalized fedvr basis
!
  use, intrinsic :: iso_c_binding
  use mod_sph, only : lmax1
  use mod_ormas, only : nfun
  use mod_const, only : czero
  use mod_hprod, only : fmat
  use mod_rad, only : nrad, nradfc

  implicit none
  logical(c_bool), intent(in) :: inner
  complex(c_double_complex), intent(in) :: wfn(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(inout) :: hwfn(1:(nrad-1), 0:lmax1, 1:nfun)

  integer(c_int), external :: util_omp_nproc
  integer(c_int), external :: util_omp_iproc
  integer(c_int) :: nproc, iproc, llr, ulr
  complex(c_double_complex), allocatable :: fmatp(:,:,:)

  nproc = util_omp_nproc()
  allocate(fmatp(1:nfun, 1:nfun, 0:(nproc-1)))

  !$omp parallel default(shared) private(iproc, llr, ulr)
  !###########################
  iproc = util_omp_iproc()
  fmatp(1:nfun, 1:nfun, iproc) = czero
  call util_omp_disp(1, nradfc, llr, ulr)
  call hprod_projallp_fmat(wfn, hwfn, fmatp(1,1,iproc), llr, ulr)
  !###########################
  !$omp end parallel

  call util_omp_gather(nfun * nfun, fmatp, fmat)
  if (inner) hwfn(1:(nrad-1), 0:lmax1, 1:nfun) = czero

  !$omp parallel default(shared) private(llr, ulr)
  !###########################
! call util_omp_disp(1, nrad - 1, llr, ulr)
  call util_omp_disp(1, nradfc, llr, ulr)
  call hprod_projallp_proj(inner, wfn, fmat, hwfn, llr, ulr)
  !###########################
  !$omp end parallel

!  fmat(1:nfun, 1:nfun) = czero
!  do iproc = 0, nproc - 1
!     do ifun = 1, nfun
!        do jfun = 1, nfun
!           fmat (jfun, ifun) = &
!         & fmat (jfun, ifun) + &
!         & fmatp(jfun, ifun, iproc)
!        end do
!     end do
!  end do

  !DEBUG
  ! write(6, "('hprod_projall: fmat (R)')")
  ! do ifun = 1, nfun
  !    do jfun = 1, nfun
  !       write(6, "(f10.5)", advance = 'no') dble(fmat(jfun, ifun))
  !    end do
  !    write(6, *)
  ! end do
  ! write(6, "('hprod_projall: fmat (I)')")
  ! do ifun = 1, nfun
  !    do jfun = 1, nfun
  !       write(6, "(f10.5)", advance = 'no') aimag(fmat(jfun, ifun))
  !    end do
  !    write(6, *)
  ! end do
  !DEBUG

  deallocate(fmatp)

end subroutine hprod_projall
!######################################################################
subroutine hprod_projallp_fmat(wfn, hwfn, fmat, llr, ulr)

  use, intrinsic :: iso_c_binding
  use mod_sph, only : lmax1
  use mod_bas, only : mval
  use mod_rad, only : nrad
  use mod_const, only : czero
  use mod_ormas, only : nfun

  implicit none
  integer(c_int), intent(in) :: llr, ulr
  complex(c_double_complex), intent(in) :: wfn(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: hwfn(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(inout) :: fmat(1:nfun, 1:nfun)

  complex(c_double_complex) :: tmp
  integer(c_int) :: ifun, jfun, l, irad

  do ifun = 1, nfun
     do jfun = 1, nfun
        if (mval(jfun) == mval(ifun)) then
           tmp = czero
           do l = 0, lmax1
              do irad = llr, ulr
                 tmp = tmp + conjg(wfn(irad, l, ifun)) * hwfn(irad, l, jfun)
              end do
           end do
           fmat(ifun, jfun) = fmat(ifun, jfun) + tmp
        end if
     end do
  end do

end subroutine hprod_projallp_fmat
!######################################################################
subroutine hprod_projallp_proj(inner, wfn, fmat, hwfn, llr, ulr)

  use, intrinsic :: iso_c_binding
  use mod_sph, only : lmax1
  use mod_bas, only : mval
  use mod_const, only : czero
  use mod_rad, only : nrad, wrad
  use mod_ormas, only : nfun, nfcore

  implicit none
  logical(c_bool), intent(in) :: inner
  integer(c_int), intent(in) :: llr, ulr
  complex(c_double_complex), intent(in) :: wfn(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: fmat(1:nfun, 1:nfun)
  complex(c_double_complex), intent(inout) :: hwfn(1:(nrad-1), 0:lmax1, 1:nfun)
  integer(c_int) :: l, irad, ifun, jfun

  if (.not. inner) then
     do ifun = nfcore + 1, nfun
        do jfun = 1, nfun
           if (mval(jfun) == mval(ifun)) then
              do l = 0, lmax1
                 do irad = llr, ulr
                    hwfn(irad, l, ifun) = hwfn(irad, l, ifun) - wfn(irad, l, jfun) * fmat(jfun, ifun)
                 end do
              end do
           end if
        end do
     end do
  else
     do ifun = nfcore + 1, nfun
        do jfun = 1, nfun
           if (mval(jfun) == mval(ifun)) then
              do l = 0, lmax1
                 do irad = llr, ulr
                    hwfn(irad, l, ifun) = hwfn(irad, l, ifun) + wfn(irad, l, jfun) * fmat(jfun, ifun)
                 end do
              end do
           end if
        end do
     end do
  end if

end subroutine hprod_projallp_proj
!######################################################################
