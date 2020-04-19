!######################################################################
subroutine hprod_mkint2_x3j(rho2, v2sph)

  use, intrinsic :: iso_c_binding
  use mod_const, only : czero
  use mod_bas, only : nbas2, mval
  use mod_hprod, only : int2e
  use mod_ormas, only : nact, ncore, nfun

  implicit none
  complex(c_double_complex), intent(in) :: rho2 (1:nbas2, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(in) :: v2sph(1:nbas2, 1:nfun, 1:nfun)

  integer(c_long), external :: util_omp_nproc
  integer(c_long), external :: util_omp_iproc
  integer(c_long) :: nproc, iproc, ifun, jfun, kfun, lfun, iact, jact, kact, lact, ibas, llb, ulb
  complex(c_double_complex) :: tmp
  complex(c_double_complex), allocatable :: int2p(:,:,:,:,:)
  complex(c_double_complex), external :: zdotu

  int2e(1:nact, 1:nact, 1:nact, 1:nact) = czero

  nproc = util_omp_nproc()
  allocate(int2p(1:nact, 1:nact, 1:nact, 1:nact, 0:(nproc-1)))

  !$omp parallel default(shared) private(iproc, ifun, jfun, kfun, lfun, llb, ulb, tmp)
  !###########################
  iproc = util_omp_iproc()
  int2p(1:nact, 1:nact, 1:nact, 1:nact, iproc) = czero
  call util_omp_disp(1, nbas2, llb, ulb)
  do iact = 1, nact
     ifun = ncore + iact
     do jact = 1, iact
        jfun = ncore + jact
        do kact = 1, nact
           kfun = ncore + kact
           do lact = 1, nact
              lfun = ncore + lact
              if (mval(ifun) + mval(kfun) == mval(jfun) + mval(lfun)) then
!                 int2p(iact,jact,kact,lact,iproc) = zdotu(ulb-llb+1,rho2(llb,ifun,jfun),1,v2sph(llb,kfun,lfun),1)

                 tmp = zdotu(ulb-llb+1,rho2(llb,ifun,jfun),1,v2sph(llb,kfun,lfun),1)
                 int2p(iact,jact,kact,lact,iproc) = tmp*(-1)**(mval(ifun)-mval(jfun))


!def                 tmp = czero
!def                 do ibas = llb, ulb
!def                    tmp = tmp + rho2(ibas, ifun, jfun) * v2sph(ibas, kfun, lfun)
!def                 end do
!def                 int2p(iact, jact, kact, lact, iproc) = tmp
              end if
           end do
        end do
     end do
  end do
  !###########################
  !$omp end parallel

  do iproc = 0, nproc - 1
     do iact = 1, nact
        do jact = 1, iact
           do kact = 1, nact
              do lact = 1, nact
                 int2e(iact, jact, kact, lact) = &
                 int2e(iact, jact, kact, lact) + &
                 int2p(iact, jact, kact, lact, iproc)
              end do
           end do
        end do
     end do
  end do
  do iact = 1, nact
     do jact = iact + 1, nact
        do kact = 1, nact
           do lact = 1, nact
              int2e(iact, jact, kact, lact) = conjg(int2e(jact, iact, lact, kact))
           end do
        end do
     end do
  end do

  deallocate(int2p)

!DEBUG
!  write(6, "('int2e:')")
!  do iact = 1, nact
!     do jact = 1, nact
!        do kact = 1, nact
!           do lact = 1, nact
!              write(6, "(4i5, 2f20.10)") iact, jact, kact, lact, int2e(iact, jact, kact, lact)
!           end do
!        end do
!     end do
!  end do
!DEBUG

end subroutine hprod_mkint2_x3j
!######################################################################
