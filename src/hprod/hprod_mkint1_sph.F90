!######################################################################
subroutine hprod_mkint1_sph(doall, orb, h0orb, h1orb, gorb)

  use, intrinsic :: iso_c_binding
  use mod_const, only : czero
  use mod_bas, only : nbas, mval
  use mod_control, only : ioorot
  use mod_hprod, only : int1e
  use mod_ormas, only : nact, ncore, nfun
  use mod_rad, only : ecs_flag

  implicit none
  logical(c_bool), intent(in) :: doall
  complex(c_double_complex), intent(in) ::   orb(1:nbas, 1:nfun)
  complex(c_double_complex), intent(in) :: h0orb(1:nbas, 1:nfun)
  complex(c_double_complex), intent(in) :: h1orb(1:nbas, 1:nfun)
  complex(c_double_complex), intent(in) ::  gorb(1:nbas, 1:nfun)

  integer(c_long), external :: util_omp_nproc
  integer(c_long), external :: util_omp_iproc
  complex(c_double_complex) :: tmp, h01
  complex(c_double_complex), allocatable :: int1p(:,:,:)
  integer(c_long) :: nproc, iproc, iact, jact, ifun, jfun, ibas, llb, ulb

! Orimo_ECS
  if (ecs_flag == 1) then
! Sato_ECS
!    call hprod_mkint1_sph_ecs(doall, orb, h0orb, h1orb, gorb)
     call hprod_mkint1_sph_ecs2(doall, orb, h0orb, h1orb, gorb)
! Sato_ECS
     return
  end if
! Orimo_ECS

  int1e(1:nact, 1:nact) = czero

  nproc = util_omp_nproc()
  allocate(int1p(1:nact, 1:nact, 0:(nproc-1)))

  !$omp parallel default(shared) private(iproc, ifun, jfun, llb, ulb, tmp, h01)
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
           if (doall .or. ioorot == 0) then
              ! <j|h0+h1+h2|i>
              do ibas = llb, ulb
                 h01 = gorb(ibas, ifun) + h1orb(ibas, ifun) + h0orb(ibas, ifun)
                 tmp = tmp + conjg(orb(ibas, jfun)) * h01
              end do
           else if (ioorot == 2) then
              ! <j|h1+h2|i>
              do ibas = llb, ulb
                 h01 = gorb(ibas, ifun) + h1orb(ibas, ifun)
                 tmp = tmp + conjg(orb(ibas, jfun)) * h01
              end do
           else if (ioorot == 1) then
              ! <j|h2|i>
              do ibas = llb, ulb
                 h01 = gorb(ibas, ifun)
                 tmp = tmp + conjg(orb(ibas, jfun)) * h01
              end do
           end if
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

!DEBUG
!  write(6, "('int1e:')")
!  do iact = 1, nact
!     do jact = 1, nact
!        write(6, "(2i5, 2f20.10)") iact, jact, int1e(iact, jact)
!     end do
!  end do
!DEBUG

end subroutine hprod_mkint1_sph
!######################################################################
!######################################################################
subroutine hprod_mkint1_sph_inner(rion, orb, h0orb, h1orb, gorb)

  use, intrinsic :: iso_c_binding
  use mod_const, only : czero
  use mod_bas, only : mval
  use mod_sph, only : lmax1
  use mod_control, only : ioorot
  use mod_hprod, only : int1e
  use mod_ormas, only : nact, ncore, nfun
  use mod_rad, only  : bra_wrad, nrad, xrad

  implicit none
  real(c_double), intent(in) :: rion
  complex(c_double_complex), intent(in) ::   orb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: h0orb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: h1orb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) ::  gorb(1:(nrad-1), 0:lmax1, 1:nfun)

  integer(c_long), external :: util_omp_nproc
  integer(c_long), external :: util_omp_iproc
  complex(c_double_complex) :: tmp, h01
  complex(c_double_complex), allocatable :: int1p(:,:,:)
  integer(c_long) :: nproc, iproc, iact, jact, ifun, jfun, l, irad, max_irad, llr, ulr


  int1e(1:nact, 1:nact) = czero

  nproc = util_omp_nproc()
  allocate(int1p(1:nact, 1:nact, 0:(nproc-1)))

  ! determine max_irad
  max_irad = nrad - 1
  do irad = 1, nrad - 1
     if (xrad(irad) > rion) then
        max_irad = irad
        exit
     end if
  end do

  !$omp parallel default(shared) private(iproc, ifun, jfun, llr, ulr, tmp, h01)
  !###########################
  iproc = util_omp_iproc()
  int1p(1:nact, 1:nact, iproc) = czero
!  call util_omp_disp(1, nrad-1, llr, ulr)
  call util_omp_disp(1, max_irad, llr, ulr)

  do iact = 1, nact
     ifun = ncore + iact
!not hermitial     do jact = 1, iact
     do jact = 1, nact
        jfun = ncore + jact

        tmp = czero
        if (mval(jfun) == mval(ifun)) then
           if (ioorot == 0) then
              ! <j|h0+h1+h2|i>
              do l = abs(mval(ifun)), lmax1
                 do irad = llr, ulr
                    h01 = gorb(irad, l, ifun) + h1orb(irad, l, ifun) + h0orb(irad, l, ifun)
                    tmp = tmp + conjg(orb(irad, l, jfun))* bra_wrad(irad) * h01
                 end do
              end do
           else if (ioorot == 2) then
              ! <j|h1+h2|i>
              do l = abs(mval(ifun)), lmax1
                 do irad = llr, ulr
                    h01 = gorb(irad, l, ifun) + h1orb(irad, l, ifun)
                    tmp = tmp + conjg(orb(irad, l, jfun)) * bra_wrad(irad) * h01
                 end do
              end do
           else if (ioorot == 1) then
              ! <j|h2|i>
              do l = abs(mval(ifun)), lmax1
                 do irad = llr, ulr
                    h01 = gorb(irad, l, ifun) 
                    tmp = tmp + conjg(orb(irad, l, jfun)) * bra_wrad(irad) * h01
                 end do
              end do
           end if
        end if
        int1p(jact, iact, iproc) = tmp
     end do
  end do
  !###########################
  !$omp end parallel

  do iproc = 0, nproc - 1
     do iact = 1, nact
        do jact = 1, nact
           int1e(jact, iact) = int1e(jact, iact) + int1p(jact, iact, iproc)
        end do
     end do
  end do

!not hermitian  do iproc = 0, nproc - 1
!not hermitian     do iact = 1, nact
!not hermitian        do jact = 1, iact
!not hermitian           int1e(jact, iact) = int1e(jact, iact) + int1p(jact, iact, iproc)
!not hermitian        end do
!not hermitian     end do
!not hermitian  end do
!not hermitian  do iact = 1, nact
!not hermitian     do jact = iact + 1, nact
!not hermitian        int1e(jact, iact) = conjg(int1e(iact, jact))
!not hermitian     end do
!not hermitian  end do

  deallocate(int1p)

end subroutine hprod_mkint1_sph_inner
!######################################################################
! Orimo_ECS
subroutine hprod_mkint1_sph_ecs(doall, orb, h0orb, h1orb, gorb)

  use, intrinsic :: iso_c_binding
  use mod_const, only : czero
  use mod_bas, only : mval
  use mod_sph, only : lmax1
  use mod_control, only : ioorot
  use mod_hprod, only : int1e
  use mod_ormas, only : nact, ncore, nfun
  use mod_rad, only  : irad_ecs, nrad

  implicit none
  logical(c_bool), intent(in) :: doall
  complex(c_double_complex), intent(in) ::   orb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: h0orb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: h1orb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) ::  gorb(1:(nrad-1), 0:lmax1, 1:nfun)

  integer(c_long), external :: util_omp_nproc
  integer(c_long), external :: util_omp_iproc
  complex(c_double_complex) :: tmp, h01
  complex(c_double_complex), allocatable :: int1p(:,:,:)
  integer(c_long) :: nproc, iproc, iact, jact, ifun, jfun, l, irad, llr, ulr


  int1e(1:nact, 1:nact) = czero

  nproc = util_omp_nproc()
  allocate(int1p(1:nact, 1:nact, 0:(nproc-1)))

  !$omp parallel default(shared) private(iproc, ifun, jfun, llr, ulr, tmp, h01)
  !###########################
  iproc = util_omp_iproc()
  int1p(1:nact, 1:nact, iproc) = czero
  call util_omp_disp(1, irad_ecs - 1, llr, ulr)

  do iact = 1, nact
     ifun = ncore + iact
!not hermitial     do jact = 1, iact
     do jact = 1, nact
        jfun = ncore + jact

        tmp = czero
        if (mval(jfun) == mval(ifun)) then
           if (doall .or. ioorot == 0) then
              ! <j|h0+h1+h2|i>
              do l = abs(mval(ifun)), lmax1
                 do irad = llr, ulr
                    h01 = gorb(irad, l, ifun) + h1orb(irad, l, ifun) + h0orb(irad, l, ifun)
                    tmp = tmp + conjg(orb(irad, l, jfun))* h01
                 end do
              end do
           else if (ioorot == 2) then
              ! <j|h1+h2|i>
              do l = abs(mval(ifun)), lmax1
                 do irad = llr, ulr
                    h01 = gorb(irad, l, ifun) + h1orb(irad, l, ifun)
                    tmp = tmp + conjg(orb(irad, l, jfun)) *  h01
                 end do
              end do
           else if (ioorot == 1) then
              ! <j|h2|i>
              do l = abs(mval(ifun)), lmax1
                 do irad = llr, ulr
                    h01 = gorb(irad, l, ifun) 
                    tmp = tmp + conjg(orb(irad, l, jfun)) * h01
                 end do
              end do
           end if
        end if
        int1p(jact, iact, iproc) = tmp
     end do
  end do
  !###########################
  !$omp end parallel

  do iproc = 0, nproc - 1
     do iact = 1, nact
        do jact = 1, nact
           int1e(jact, iact) = int1e(jact, iact) + int1p(jact, iact, iproc)
        end do
     end do
  end do
!!not hermitian  do iproc = 0, nproc - 1
!!not hermitian     do iact = 1, nact
!!not hermitian        do jact = 1, iact
!!not hermitian           int1e(jact, iact) = int1e(jact, iact) + int1p(jact, iact, iproc)
!!not hermitian        end do
!!not hermitian     end do
!!not hermitian  end do
!!not hermitian  do iact = 1, nact
!!not hermitian     do jact = iact + 1, nact
!!not hermitian        int1e(jact, iact) = conjg(int1e(iact, jact))
!!not hermitian     end do
!!not hermitian  end do

  deallocate(int1p)

end subroutine hprod_mkint1_sph_ecs
! Orimo_ECS
!######################################################################
subroutine hprod_mkint1_sph_ecs2(doall, orb, h0orb, h1orb, gorb)

  use, intrinsic :: iso_c_binding
  use mod_const, only : czero
  use mod_bas, only : mval, tmat
  use mod_control, only : ioorot
  use mod_hprod, only : int1e, h0orb_out, h1orb_out
  use mod_ormas, only : nact, ncore, nfun
  use mod_rad, only : irad_ecs, nrad, ndvr
  use mod_sph, only : lmax1

  implicit none
  logical(c_bool), intent(in) :: doall
  complex(c_double_complex), intent(in) ::   orb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: h0orb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: h1orb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) ::  gorb(1:(nrad-1), 0:lmax1, 1:nfun)

  integer(c_long), external :: util_omp_nproc
  integer(c_long), external :: util_omp_iproc
  complex(c_double_complex) :: tmp, h01
  complex(c_double_complex), allocatable :: int1p(:,:,:)
  integer(c_long) :: nproc, iproc, iact, jact, ifun, jfun, irad, l, llr, ulr, jrad

  int1e(1:nact, 1:nact) = czero
  nproc = util_omp_nproc()
  allocate(int1p(1:nact, 1:nact, 0:(nproc-1)))

  !$omp parallel default(shared) private(iproc, ifun, jfun, llr, ulr, tmp, h01, jrad)
  !###########################
  iproc = util_omp_iproc()
  int1p(1:nact, 1:nact, iproc) = czero
!  call util_omp_disp(1, nrad-1, llr, ulr)
  call util_omp_disp(1, irad_ecs-1, llr, ulr)

  do iact = 1, nact
     ifun = ncore + iact
     do jact = 1, iact
        jfun = ncore + jact

        tmp = czero
        if (mval(jfun) == mval(ifun)) then
           if (doall .or. ioorot == 0) then
              ! <j|h0+h1+h2|i>
              do l = abs(mval(ifun)), lmax1
                 do irad = llr, ulr
                    h01 = gorb(irad, l, ifun) + h1orb(irad, l, ifun) + h0orb(irad, l, ifun)
                    tmp = tmp + conjg(orb(irad, l, jfun)) * h01
! Sato_ECS
                    h01 = h1orb_out(irad,l,jfun) + h0orb_out(irad,l,jfun)
                    tmp = tmp + conjg(h01) * orb(irad,l,ifun)
! Sato_ECS
                 end do
              end do
           else if (ioorot == 2) then
              ! <j|h1+h2|i>
              do l = abs(mval(ifun)), lmax1
                 do irad = llr, ulr
                    h01 = gorb(irad, l, ifun) + h1orb(irad, l, ifun)
                    tmp = tmp + conjg(orb(irad, l, jfun)) * h01
! Sato_ECS
                    h01 = h1orb_out(irad,l,jfun)
                    tmp = tmp + conjg(h01) * orb(irad,l,ifun)
! Sato_ECS
                 end do
              end do
           else if (ioorot == 1) then
              ! <j|h2|i>
              do l = abs(mval(ifun)), lmax1
                 do irad = llr, ulr
                    h01 = gorb(irad, l, ifun)
                    tmp = tmp + conjg(orb(irad, l, jfun)) * h01
                 end do
              end do
           end if
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

end subroutine hprod_mkint1_sph_ecs2
!######################################################################
