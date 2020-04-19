!######################################################################
subroutine hprod_projpp(zfac, wfn, hwfn)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad
  use mod_sph, only : lmax1
  use mod_ormas, only : nfun
  use mod_control, only : psp_type
  use mod_bas, only : mval,pp_maxl,pp_nump,pp_irmax,pp_fproj,pp_gproj

  implicit none
  complex(c_double_complex), intent(in) :: zfac
  complex(c_double_complex), intent(in) :: wfn(1:(nrad-1), 0:lmax1, 1:*)
  complex(c_double_complex), intent(inout) :: hwfn(1:(nrad-1), 0:lmax1, 1:*)
  integer(c_long) :: ifun, l, i, irad
  complex(c_double_complex) :: ovlp

  if (psp_type.ne.1 .and. psp_type.ne.5) return

  !$omp parallel default(shared) private(ovlp)
  !$omp do
  do ifun = 1, nfun
     do l = abs(mval(ifun)), min(pp_maxl,lmax1)
        do i = 1, pp_nump(l)
           ovlp = 0.d+0
           do irad = 1, pp_irmax(i,l)
              ovlp = ovlp + pp_gproj(irad,i,l) * wfn(irad,l,ifun)
           end do
           ovlp = ovlp * zfac

           do irad = 1, pp_irmax(i,l)
              hwfn(irad,l,ifun) = hwfn(irad,l,ifun) + pp_fproj(irad,i,l) * ovlp
           end do
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine hprod_projpp
!######################################################################
