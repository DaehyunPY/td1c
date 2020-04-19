!#######################################################################
subroutine pes_spec1_kz(wfn, cic)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : xrad
  use mod_hprod, only : den1
  use mod_const, only : czero, ctwo
  use mod_pes, only : pes_psik, pes_rhok, pes_llr
  use mod_ormas, only : ncore, nact, nfun

  implicit none
  complex(c_double_complex), intent(in) :: wfn(*)
  complex(c_double_complex), intent(in) :: cic(*)
  complex(c_double_complex), allocatable :: dens(:,:), ovlp(:,:)
  integer(c_long) :: ifun, iact, jact

  allocate(ovlp(1:nfun, 1:nfun))
  allocate(dens(1:nfun, 1:nfun))
  ovlp(1:nfun, 1:nfun) = czero
  dens(1:nfun, 1:nfun) = czero

  if (.true.) then
     call ormas_mkden1(cic, den1)
     do ifun = 1, ncore
        dens(ifun, ifun) = ctwo
     end do
     do iact = 1, nact
!DEBUG        dens(ncore+iact, ncore+iact) = den1(iact,iact)
        do jact = 1, nact
           dens(ncore + jact, ncore + iact) = den1(jact, iact)
        end do
     end do
  else
!!DEBUG
!     call ormas_mkden1(cic, den1)
!!DEBUG
!     call hprod_mkovlp(xrad(pes_llr), wfn, wfn, ovlp)
!     call ormas_mkidm(ovlp, cic, dens)
  end if

  call pes_spec1_k_psik(wfn, pes_psik)
  call pes_spec1_k_rhok_z(pes_psik, pes_rhok)
  call pes_spec1_k_trace(dens, pes_rhok)

  deallocate(dens)
  deallocate(ovlp)

end subroutine pes_spec1_kz
!#######################################################################
subroutine pes_spec1_k_rhok_z(psik, rhok)

  use, intrinsic :: iso_c_binding
  use mod_const, only : czero, pi, iunit, one, two
  use mod_pes, only : pes_numk
  use mod_sph, only : lmax1
  use mod_ormas, only : nfcore, nfun
  use mod_bas, only : mval

  implicit none
  complex(c_double_complex), intent(in) :: psik(0:pes_numk, 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(out) :: rhok(0:pes_numk, 1:nfun, 1:nfun)

  integer(c_long) :: ifun, jfun, ik, l
  integer(c_long) :: iproc, llk, ulk
  integer(c_long), external :: util_omp_iproc
  complex(c_double_complex) :: psizk(0:pes_numk, 1:nfun)
  complex(c_double_complex) :: temp

  rhok(0:pes_numk, 1:nfun, 1:nfun) = czero
  psizk(0:pes_numk, 1:nfun) = czero
  temp = two * sqrt(pi)

  !$omp parallel default(shared) private(iproc, llk, ulk)
  !###########################
  iproc = util_omp_iproc()
  call util_omp_disp(0, pes_numk, llk, ulk)
  do ifun = nfcore + 1, nfun
     do ik = llk, ulk
        do l = 0, lmax1
          if (mval(ifun) == 0) then
             psizk(ik, ifun) = psizk(ik, ifun) + ((- iunit) ** l) * sqrt((two * l + one) * one) * psik(ik, l, ifun) 
          end if
        end do
        psizk(ik, ifun) = psizk(ik, ifun) * temp
     end do
  end do
  do ifun = nfcore + 1, nfun
     do jfun = nfcore + 1, nfun
        do ik = llk, ulk
           rhok(ik, ifun, jfun) = conjg(psizk(ik, ifun)) * psizk(ik, jfun)
        end do
     end do
  end do
  !###########################
  !$omp end parallel
  
end subroutine pes_spec1_k_rhok_z
!#######################################################################
