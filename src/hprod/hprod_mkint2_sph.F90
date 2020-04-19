!######################################################################
subroutine hprod_mkint2_sph(rho2, v2sph)

  use, intrinsic :: iso_c_binding
  use mod_const, only : czero
  use mod_bas, only : nbas2, mval
  use mod_hprod, only : int2e
  use mod_ormas, only : nact, ncore, nfun
  use mod_rad, only : ecs_flag

  implicit none
  complex(c_double_complex), intent(in) :: rho2 (1:nbas2, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(in) :: v2sph(1:nbas2, 1:nfun, 1:nfun)

  integer(c_long), external :: util_omp_nproc
  integer(c_long), external :: util_omp_iproc
  integer(c_long) :: nproc, iproc, ifun, jfun, kfun, lfun, iact, jact, kact, lact, ibas, llb, ulb
  complex(c_double_complex) :: tmp
  complex(c_double_complex), allocatable :: int2p(:,:,:,:,:)
  complex(c_double_complex), external :: zdotu

  if (ecs_flag == 1) then
! Sato_ECS
!    call hprod_mkint2_sph_ecs(rho2, v2sph)
     call hprod_mkint2_sph_ecs2(rho2, v2sph)
! Sato_ECS
     return
  end if

  nproc = util_omp_nproc()
  int2e(1:nact, 1:nact, 1:nact, 1:nact) = czero
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
                 int2p(iact,jact,kact,lact,iproc) = zdotu(ulb-llb+1,rho2(llb,ifun,jfun),1,v2sph(llb,kfun,lfun),1)
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

end subroutine hprod_mkint2_sph
!######################################################################
subroutine hprod_mkint2_sph_ecs(rho2, v2sph)

  use, intrinsic :: iso_c_binding
  use mod_const, only : czero
  use mod_bas, only : nbas2, mval
  use mod_hprod, only : int2e
  use mod_ormas, only : nact, ncore, nfun
  use mod_rad, only : irad_ecs, nrad
  use mod_sph, only : lmax2

  implicit none
  complex(c_double_complex), intent(in) :: rho2 (1:(nrad-1), 0:lmax2, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(in) :: v2sph(1:(nrad-1), 0:lmax2, 1:nfun, 1:nfun)

  integer(c_long), external :: util_omp_nproc
  integer(c_long), external :: util_omp_iproc
  integer(c_long) :: nproc, iproc, ifun, jfun, kfun, lfun, iact, jact, kact, lact, ibas, llb, ulb
  complex(c_double_complex) :: tmp
  complex(c_double_complex), allocatable :: int2p(:,:,:,:,:)
  complex(c_double_complex), external :: zdotu

  integer(c_long) :: irad, l, llr, ulr


  int2e(1:nact, 1:nact, 1:nact, 1:nact) = czero

  nproc = util_omp_nproc()
  allocate(int2p(1:nact, 1:nact, 1:nact, 1:nact, 0:(nproc-1)))

  !$omp parallel default(shared) private(iproc, ifun, jfun, kfun, lfun, llr, ulr, tmp)
  !###########################
  iproc = util_omp_iproc()
  int2p(1:nact, 1:nact, 1:nact, 1:nact, iproc) = czero
  call util_omp_disp(1, irad_ecs-1, llr, ulr) !! <<<<<<< changed
  do iact = 1, nact
     ifun = ncore + iact
     do jact = 1, nact !! <<<<<<< changed! Sato_ECS
        jfun = ncore + jact
        do kact = 1, nact
           kfun = ncore + kact
           do lact = 1, nact
              lfun = ncore + lact
              if (mval(ifun) + mval(kfun) == mval(jfun) + mval(lfun)) then
                 !def int2p(iact,jact,kact,lact,iproc) = zdotu(ulb-llb+1,rho2(llb,ifun,jfun),1,v2sph(llb,kfun,lfun),1)
                 tmp = czero
                 do l = 0, lmax2
                    do irad = llr, ulr
                       tmp = tmp + rho2(irad, l , ifun, jfun) * v2sph(irad, l, kfun, lfun)
                    end do
                 end do
                 int2p(iact, jact, kact, lact, iproc) = tmp
                 
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
        do jact = 1, nact !! <<<<<< changed
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
!not hermitian  do iact = 1, nact
!not hermitian     do jact = iact + 1, nact
!not hermitian        do kact = 1, nact
!not hermitian           do lact = 1, nact
!not hermitian              int2e(iact, jact, kact, lact) = conjg(int2e(jact, iact, lact, kact))
!not hermitian           end do
!not hermitian        end do
!not hermitian     end do
!not hermitian  end do

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

end subroutine hprod_mkint2_sph_ecs
!######################################################################
subroutine hprod_mkint2_sph_ecs2(rho2, v2sph)

  use, intrinsic :: iso_c_binding
  use mod_const, only : czero
  use mod_bas, only : nbas2, mval
  use mod_hprod, only : int2e
  use mod_ormas, only : nact, ncore, nfun
  use mod_rad, only : irad_ecs, nrad
  use mod_sph, only : lmax2

  implicit none
  complex(c_double_complex), intent(in) :: rho2 (1:(nrad-1), 0:lmax2, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(in) :: v2sph(1:(nrad-1), 0:lmax2, 1:nfun, 1:nfun)

  integer(c_long), external :: util_omp_nproc
  integer(c_long), external :: util_omp_iproc
  integer(c_long) :: nproc, iproc, ifun, jfun, kfun, lfun, iact, jact, kact, lact, ibas, llb, ulb
  complex(c_double_complex) :: tmp
  complex(c_double_complex), allocatable :: int2p(:,:,:,:,:)
  complex(c_double_complex), external :: zdotu

  integer(c_long) :: irad, l, llr, ulr


  int2e(1:nact, 1:nact, 1:nact, 1:nact) = czero

  nproc = util_omp_nproc()
  allocate(int2p(1:nact, 1:nact, 1:nact, 1:nact, 0:(nproc-1)))

  !$omp parallel default(shared) private(iproc, ifun, jfun, kfun, lfun, llr, ulr, tmp)
  !###########################
  iproc = util_omp_iproc()
  int2p(1:nact, 1:nact, 1:nact, 1:nact, iproc) = czero
  call util_omp_disp(1, irad_ecs-1, llr, ulr) !! <<<<<<< changed
  do iact = 1, nact
     ifun = ncore + iact
     do jact = 1, iact
        jfun = ncore + jact
        do kact = 1, nact
           kfun = ncore + kact
           do lact = 1, nact
              lfun = ncore + lact
              if (mval(ifun) + mval(kfun) == mval(jfun) + mval(lfun)) then
                 tmp = czero
                 do l = 0, lmax2
                    do irad = llr, ulr
                       tmp = tmp + rho2(irad, l , ifun, jfun) * v2sph(irad, l, kfun, lfun)
                    end do
                 end do
                 int2p(iact, jact, kact, lact, iproc) = tmp
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

end subroutine hprod_mkint2_sph_ecs2
!######################################################################
