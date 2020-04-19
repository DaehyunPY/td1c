!######################################################################
subroutine hprod_projpp(zfac, lfield, wfn, hwfn)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad
  use mod_sph, only : lmax1, nlat
  use mod_ormas, only : nfun
  use mod_control, only : icomp, psp_type, igauge
  use mod_bas, only : mval,pp_maxl,pp_nump,pp_irmax,pp_fproj,pp_gproj

  implicit none
  complex(c_double_complex), intent(in) :: zfac
  real(c_double), intent(in) :: lfield(1:3, 1:3)
  complex(c_double_complex), intent(in) :: wfn(1:(nrad-1), 0:lmax1, 1:*)
  complex(c_double_complex), intent(inout) :: hwfn(1:(nrad-1), 0:lmax1, 1:*)
  integer(c_int) :: ifun, l, i, irad, max_irmax
  complex(c_double_complex) :: ovlp
  complex(c_double_complex), allocatable :: winp(:,:,:)
  complex(c_double_complex), allocatable :: wout(:,:,:)

  if (psp_type.ne.1 .and. psp_type.ne.5 .and. psp_type.ne.6 .and. psp_type.ne.7) return

!v5.7  !$omp parallel default(shared) private(ovlp)
!v5.7  !$omp do
!v5.7  do ifun = 1, nfun
!v5.7     do l = abs(mval(ifun)), min(pp_maxl,lmax1)
!v5.7        do i = 1, pp_nump(l)
!v5.7           ovlp = 0.d+0
!v5.7           do irad = 1, pp_irmax(i,l)
!v5.7              ovlp = ovlp + pp_gproj(irad,i,l) * wfn(irad,l,ifun)
!v5.7           end do
!v5.7           ovlp = ovlp * zfac
!v5.7
!v5.7           do irad = 1, pp_irmax(i,l)
!v5.7              hwfn(irad,l,ifun) = hwfn(irad,l,ifun) + pp_fproj(irad,i,l) * ovlp
!v5.7           end do
!v5.7        end do
!v5.7     end do
!v5.7  end do
!v5.7  !$omp end do
!v5.7  !$omp end parallel

  max_irmax = 0
  do l = 0, min(pp_maxl,lmax1)
     max_irmax = max(max_irmax,maxval(pp_irmax(1:pp_nump(l),l)))
  end do

  allocate(winp(1:max_irmax,0:lmax1,1:nfun))
  allocate(wout(1:max_irmax,0:lmax1,1:nfun))

  winp = 0D0
  wout = 0D0
  if (icomp == 1 .and. igauge == 1) then
     do ifun = 1, nfun
        !do l = abs(mval(ifun)), min(pp_maxl,lmax1)
        do l = abs(mval(ifun)), lmax1
           do irad = 1, max_irmax
              wout(irad,l,ifun) = wfn(irad,l,ifun)
           end do
        end do
     end do
     call hprod_gtran(.false., lfield, max_irmax, wout, winp)
  else
     do ifun = 1, nfun
        do l = abs(mval(ifun)), min(pp_maxl,lmax1)
           do irad = 1, max_irmax
              winp(irad,l,ifun) = wfn(irad,l,ifun)
           end do
        end do
     end do
  end if

  wout = 0D0
  !$omp parallel default(shared) private(ovlp)
  !$omp do
  do ifun = 1, nfun
     do l = abs(mval(ifun)), min(pp_maxl,lmax1)
        do i = 1, pp_nump(l)
           ovlp = 0.d+0
           do irad = 1, pp_irmax(i,l)
              !ovlp = ovlp + pp_gproj(irad,i,l) * wfn(irad,l,ifun)
              ovlp = ovlp + pp_gproj(irad,i,l) * winp(irad,l,ifun)
           end do
           ovlp = ovlp * zfac

           do irad = 1, pp_irmax(i,l)
              !hwfn(irad,l,ifun) = hwfn(irad,l,ifun) + pp_fproj(irad,i,l) * ovlp
              wout(irad,l,ifun) = wout(irad,l,ifun) + pp_fproj(irad,i,l) * ovlp
           end do
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel

  if (icomp == 1 .and. igauge == 1) then
     call hprod_gtran(.true., lfield, max_irmax, wout, winp)
     do ifun = 1, nfun
        !do l = abs(mval(ifun)), min(pp_maxl,lmax1)
        do l = abs(mval(ifun)), lmax1
           do irad = 1, max_irmax
              hwfn(irad,l,ifun) = hwfn(irad,l,ifun) + winp(irad,l,ifun)
           end do
        end do
     end do
  else
     do ifun = 1, nfun
        do l = abs(mval(ifun)), min(pp_maxl,lmax1)
           do irad = 1, max_irmax
              hwfn(irad,l,ifun) = hwfn(irad,l,ifun) + wout(irad,l,ifun)
           end do
        end do
     end do
  end if

  deallocate(wout)
  deallocate(winp)

end subroutine hprod_projpp
!######################################################################
