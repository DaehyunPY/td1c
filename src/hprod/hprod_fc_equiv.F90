!######################################################################
subroutine hprod_fc_equiv(wfn)

  use, intrinsic :: iso_c_binding
  use mod_sph, only : lmax1
  use mod_ormas, only : nfun, nfcore
  use mod_bas, only : nval, lval, mval
  use mod_rad, only : nrad
  use mod_const, only : czero

  implicit none
  complex(c_double_complex), intent(inout) :: wfn(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), allocatable :: twfn(:,:)
  integer(c_long) :: llr, ulr, ifun, jfun

  if (nfcore == 0) return
  write(6, "('# FC orbitals are explicitly spherical-symmetry adapted.')")
  allocate(twfn(1:(nrad-1), 1:nfcore))

  !$omp parallel default(shared) private(llr, ulr)
  !###########################
  call util_omp_disp(1, nrad-1, llr, ulr)
  do ifun = 1, nfcore
     twfn(llr:ulr, ifun) = wfn(llr:ulr, lval(ifun), ifun)
  end do
  wfn(llr:ulr, 0:lmax1, 1:nfcore) = czero

  do ifun = 1, nfcore
     if (mval(ifun) == 0) then
        wfn(llr:ulr, lval(ifun), ifun) = twfn(llr:ulr, ifun)
     else
        do jfun = 1, nfcore
           if (nval(jfun) == nval(ifun) .and. lval(jfun) == lval(ifun) .and. mval(jfun) == 0) then
              wfn(llr:ulr, lval(ifun), ifun) = twfn(llr:ulr, jfun)
              exit
           end if
        end do
     end if
  end do
  !###########################
  !$omp end parallel

  deallocate(twfn)

end subroutine hprod_fc_equiv
!######################################################################
