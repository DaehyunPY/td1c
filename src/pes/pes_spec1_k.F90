!#######################################################################
subroutine pes_spec1_k(wfn, cic)

  use, intrinsic :: iso_c_binding
  use mod_hprod, only : den1
  use mod_const, only : czero, ctwo
  use mod_pes, only : pes_psik, pes_rhok
  use mod_ormas, only : ncore, nact, nfun

  implicit none
  complex(c_double_complex), intent(in) :: wfn(*)
  complex(c_double_complex), intent(in) :: cic(*)
  complex(c_double_complex), allocatable :: dens(:,:)
  integer(c_long) :: ifun, iact, jact

  allocate(dens(1:nfun, 1:nfun))
  dens(1:nfun, 1:nfun) = czero

  if (.true.) then
     call ormas_mkden1(cic, den1)
     do ifun = 1, ncore
        dens(ifun, ifun) = ctwo
     end do
     do iact = 1, nact
        do jact = 1, nact
           dens(ncore + jact, ncore + iact) = den1(jact, iact)
        end do
     end do
  end if

  call pes_spec1_k_psik(wfn, pes_psik)
  call pes_spec1_k_rhok(pes_psik, pes_rhok)
  call pes_spec1_k_trace(dens, pes_rhok)

  deallocate(dens)

end subroutine pes_spec1_k
!#######################################################################
subroutine pes_spec1_k_psik(wfn, psik)

  use, intrinsic :: iso_c_binding
  use mod_const, only : czero
  use mod_pes, only : pes_numk, pes_llr, pes_ulr, pes_k_min, pes_k_max, pes_k_step, pes_bess
  use mod_rad, only : nrad, xrad, wrad
  use mod_sph, only : lmax1
  use mod_ormas, only : nfcore, nfun

  implicit none
  complex(c_double_complex), intent(in) :: wfn(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(out) :: psik(0:pes_numk, 0:lmax1, 1:nfun)

  integer(c_long) :: ifun, irad, ik, l
  integer(c_long) :: iproc, llk, ulk
  integer(c_long), external :: util_omp_iproc
  real(c_double) :: rwval

!SAMPLE TO SHOW THE USAGE
!  write(6, "(' pes_numk = ', i5)") pes_numk
!  write(6, "(' pes_llr  = ', i5)") pes_llr
!  write(6, "(' pes_ulr  = ', i5)") pes_ulr
!  write(6, "(' pes_k_min  = ', f20.10)") pes_k_min
!  write(6, "(' pes_k_max  = ', f20.10)") pes_k_max  
!  write(6, "(' pes_k_step = ', f20.10)") pes_k_step
!  ik = 0
!  write(6, "('pes_spec1_k_psik: sph-bess for k = ', f20.10)") pes_k_min + ik * pes_k_step
!  do irad = pes_llr, pes_ulr
!     rwval = xrad(irad) * sqrt(wrad(irad))
!     write(6, "(i5, f20.10)", advance = 'no') irad, xrad(irad)
!     do l = 0, lmax1
!        write(6, "(f20.10)", advance = 'no') pes_bess(irad, ik, l) / rwval 
!     end do
!     write(6, *)
!  end do
! stop
!SAMPLE

  psik(0:pes_numk, 0:lmax1, 1:nfun) = czero

  !$omp parallel default(shared) private(iproc, llk, ulk)
  !###########################
  iproc = util_omp_iproc()
  call util_omp_disp(0, pes_numk, llk, ulk)
  do ifun = nfcore + 1, nfun
     do ik = llk, ulk
        do irad = pes_llr, pes_ulr
           do l = 0, lmax1
              psik(ik, l, ifun) = psik(ik, l, ifun) + pes_bess(irad, ik, l) * wfn(irad, l, ifun)
           end do
        end do
     end do
  end do
  !###########################
  !$omp end parallel

end subroutine pes_spec1_k_psik
!#######################################################################
subroutine pes_spec1_k_rhok(psik, rhok)
  !
  ! rhok(k, i, j) <= sum_{l} conjg(psik(k, l, i)) * psik(k, l, j)
  !
  use, intrinsic :: iso_c_binding
  use mod_const, only : czero, pi
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
  complex(c_double_complex) :: temp

  temp = 16.D0 * pi * pi
  rhok(0:pes_numk, 1:nfun, 1:nfun) = czero

  !$omp parallel default(shared) private(iproc, llk, ulk)
  !###########################
  iproc = util_omp_iproc()
  call util_omp_disp(0, pes_numk, llk, ulk)
  do ifun = nfcore + 1, nfun
     do jfun = nfcore + 1, nfun
        do ik = llk, ulk
           do l = 0, lmax1
              if (mval(ifun) == mval(jfun)) then
                 rhok(ik, ifun, jfun) = rhok(ik, ifun, jfun) + conjg(psik(ik, l, ifun)) * psik(ik, l, jfun)
              end if
           end do
           rhok(ik, ifun, jfun) = rhok(ik, ifun, jfun) * temp
        end do
     end do
  end do
  !###########################
  !$omp end parallel
  
end subroutine pes_spec1_k_rhok
!#######################################################################
subroutine pes_spec1_k_trace(dens, rhok)
  !
  ! rhok(k, 1, 1) <= sum_{ij} rhok(k, i, j) * dens(j, i)
  !
  use, intrinsic :: iso_c_binding
  use mod_const, only : czero,pi
  use mod_pes, only : pes_numk, pes_k_min, pes_k_step
  use mod_ormas, only : nfun

  implicit none
  complex(c_double_complex), intent(in) :: dens(1:nfun, 1:nfun)
  complex(c_double_complex), intent(inout) :: rhok(0:pes_numk, 1:nfun, 1:nfun)

  integer(c_long) :: ifun, jfun, irad, ik, l
  integer(c_long) :: iproc, llk, ulk
  integer(c_long), external :: util_omp_iproc
  complex(c_double_complex) :: temp
  real(c_double) :: kval 

  !$omp parallel default(shared) private(iproc, llk, ulk, kval, temp)
  !###########################
  iproc = util_omp_iproc()
  call util_omp_disp(0, pes_numk, llk, ulk)
  do ik = llk, ulk
      kval = pes_k_min + ik * pes_k_step
      temp = czero
!1RDM      do ifun = nfcore + 1, ncore
!1RDM         temp = temp + 2 * rhok(ik, ifun, ifun)
!1RDM      end do
!1RDM      do ifun = 1, nact
!1RDM         do jfun = 1, nact
!1RDM            temp = temp + rhok(ik, ncore + ifun, ncore + jfun) * dens(jfun, ifun)
!1RDM         end do
!1RDM      end do
      do ifun = 1, nfun
         do jfun = 1, nfun
            temp = temp + rhok(ik, ifun, jfun) * dens(jfun, ifun)
         end do
      end do
! Sato_tSURFF
      !rhok(ik, 1, 1) = temp
      rhok(ik, 1, 1) = temp/(2d0*pi)**3
! Sato_tSURFF
  end do
  !###########################
  !$omp end parallel

end subroutine pes_spec1_k_trace
!#######################################################################
