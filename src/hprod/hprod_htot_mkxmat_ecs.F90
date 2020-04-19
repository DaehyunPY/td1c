!######################################################################
subroutine hprod_htot_mkxmat_fc_ecs(dtime, orb, h0orb, h1orb, gorb, v2orb)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : nbas, mval
  use mod_const, only : czero, runit, iunit
  use mod_ormas, only : nfcore, nfun
  use mod_control, only : icomp
  use mod_hprod, only : xmat

  use mod_rad, only : irad_ecs, nrad ! <<<<< added
  use mod_sph, only : lmax1  ! <<<<< added

  implicit none
  real(c_double), intent(in) :: dtime
  complex(c_double_complex), intent(in) ::   orb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: h0orb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: h1orb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) ::  gorb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: v2orb(1:(nrad-1), 0:lmax1, 1:nfun)

  integer(c_int), external :: util_omp_nproc
  integer(c_int), external :: util_omp_iproc
  complex(c_double_complex) :: tmp, h01, h2, tfac
  complex(c_double_complex), allocatable :: xmatp(:,:,:)

  integer(c_int) :: nproc, iproc, ifun, jfun, l, llr, ulr, irad ! <<<<<<< changed

  if (nfcore == 0) return
  xmat(1:nfcore, (nfcore+1):nfun) = czero

  nproc = util_omp_nproc()
  allocate(xmatp(1:nfun, 1:nfun, 0:(nproc-1)))

  !$omp parallel default(shared) private(iproc, llr, ulr, tmp, h01, h2) ! <<<<<<< changed
  !###########################
  iproc = util_omp_iproc()
  xmatp(1:nfun, 1:nfun, iproc) = czero

  call util_omp_disp(1, irad_ecs-1, llr, ulr) ! <<<<<<< changed
  do ifun = nfcore + 1, nfun
     do jfun = 1, nfcore
        tmp = czero
        if (mval(jfun) == mval(ifun)) then
           ! <k|(h0+v1+zE) * |i>
           do l = abs(mval(ifun)), lmax1
              do irad = llr, ulr
                 h01 = h0orb(irad, l, jfun) + h1orb(irad, l, jfun)
                 tmp = tmp + conjg(h01) * orb(irad, l, ifun)
              end do
           end do
           !do ibas = llb, ulb
           !   h01 = h0orb(ibas, jfun) + h1orb(ibas, jfun)
           !   tmp = tmp + conjg(h01) * orb(ibas, ifun)
           !end do

           ! <k| * v2|i>
           do l = abs(mval(ifun)), lmax1
              do irad = llr, ulr
                 h2 = gorb(irad, l, ifun) + v2orb(irad, l, ifun)
                 tmp = tmp + conjg(orb(irad, l, jfun)) * h2
              end do
           end do
           !do ibas = llb, ulb
           !   h2 = gorb(ibas, ifun) + v2orb(ibas, ifun)
           !   tmp = tmp + conjg(orb(ibas, jfun)) * h2
           !end do

        end if
        xmatp(jfun, ifun, iproc) = tmp
     end do
  end do

  !###########################
  !$omp end parallel

  do iproc = 0, nproc - 1
     do ifun = nfcore + 1, nfun
        do jfun = 1, nfcore
           xmat(jfun, ifun) = xmat(jfun, ifun) + xmatp(jfun, ifun, iproc)
        end do
     end do
  end do

  if (icomp == 0) then
     tfac = runit * dtime
  else
     tfac = iunit * dtime
  end if
  xmat(1:nfcore, (nfcore+1):nfun) = xmat(1:nfcore, (nfcore+1):nfun) * tfac

  deallocate(xmatp)

end subroutine hprod_htot_mkxmat_fc_ecs
!######################################################################
subroutine hprod_htot_mkxmat_cc_ecs(dtime, orb, h0orb, h1orb, gorb, v2orb)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : nbas, mval
  use mod_const, only : czero, runit, iunit
  use mod_ormas, only : nfcore, ncore, nfun
  use mod_control, only : icomp, ioorot
  use mod_hprod, only : xmat, h0orb_out, h1orb_out

  use mod_rad, only : irad_ecs, nrad ! <<<<< added
  use mod_sph, only : lmax1  ! <<<<< added

  implicit none
  real(c_double), intent(in) :: dtime
  complex(c_double_complex), intent(in) ::   orb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: h0orb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: h1orb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) ::  gorb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: v2orb(1:(nrad-1), 0:lmax1, 1:nfun)

  integer(c_int), external :: util_omp_nproc
  integer(c_int), external :: util_omp_iproc
  complex(c_double_complex) :: tmp, h012, tfac
  complex(c_double_complex), allocatable :: xmatp(:,:,:)
  integer(c_int) :: nproc, iproc, ifun, jfun, l, irad, llr, ulr

  if (ncore == 0) return
  if (ncore == nfcore) return

  xmat((nfcore+1):ncore, (nfcore+1):ncore) = czero

  nproc = util_omp_nproc()
  allocate(xmatp(1:nfun, 1:nfun, 0:(nproc-1)))

  !$omp parallel default(shared) private(iproc, llr, ulr, tmp, h012)
  !###########################
  iproc = util_omp_iproc()
  xmatp(1:nfun, 1:nfun, iproc) = czero
  call util_omp_disp(1, irad_ecs-1, llr, ulr)

  do ifun = nfcore + 1, ncore
     do jfun = nfcore + 1, ncore
        tmp = czero
        if (mval(jfun) == mval(ifun)) then
           if (ioorot == 0) then
              ! <j|h0+h1+h2|i>
              do l = abs(mval(ifun)), lmax1
                 do irad = llr, ulr
                    h012 = v2orb(irad, l, ifun) + gorb(irad, l, ifun) + h1orb(irad, l, ifun) + h0orb(irad, l, ifun)
                    tmp = tmp + conjg(orb(irad, l, jfun)) * h012
! Sato_ECS
                    h012 = h1orb_out(irad,l,jfun) + h0orb_out(irad,l,jfun)
                    tmp = tmp + conjg(h012) * orb(irad,l,ifun)
! Sato_ECS
                 end do
              end do
           else if (ioorot == 2) then
              ! <j|h1+h2|i>
              do l = abs(mval(ifun)), lmax1
                 do irad = llr, ulr
                    h012 = v2orb(irad, l, ifun) + gorb(irad, l, ifun) + h1orb(irad, l, ifun)
                    tmp = tmp + conjg(orb(irad, l, jfun)) * h012
! Sato_ECS
                    h012 = h1orb_out(irad,l,jfun)
                    tmp = tmp + conjg(h012) * orb(irad,l,ifun)
! Sato_ECS
                 end do
              end do
           else if (ioorot == 1) then
              ! <j|h2|i>
              do l = abs(mval(ifun)), lmax1
                 do irad = llr, ulr
                    h012 = v2orb(irad, l, ifun) + gorb(irad, l, ifun)
                    tmp = tmp + conjg(orb(irad, l, jfun)) * h012
                 end do
              end do
           end if
        end if
        xmatp(jfun, ifun, iproc) = xmatp(jfun, ifun, iproc) + tmp
     end do
  end do

  !###########################
  !$omp end parallel

  do iproc = 0, nproc - 1
     do ifun = nfcore + 1, ncore
        do jfun = nfcore + 1, ncore
           xmat (jfun, ifun) = &
         & xmat (jfun, ifun) + &
         & xmatp(jfun, ifun, iproc)
        end do
     end do
  end do

  if (icomp == 0) then
     tfac = runit * dtime
  else
     tfac = iunit * dtime
  end if
  xmat((nfcore+1):ncore, (nfcore+1):ncore) = xmat((nfcore+1):ncore, (nfcore+1):ncore) * tfac

  deallocate(xmatp)

end subroutine hprod_htot_mkxmat_cc_ecs
!######################################################################
subroutine hprod_htot_mkxmat_aa1_ecs(dtime, orb, h0orb, h1orb, gorb, v2orb)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : nbas, mval, tmat
  use mod_const, only : czero, runit, iunit
  use mod_ormas, only : ncore, nact, nfun, nsub, lorb_sub
  use mod_control, only : icomp, ioorot
  use mod_hprod, only : xmat, fmat, h0orb_out, h1orb_out

  use mod_rad, only : irad_ecs, nrad, ndvr ! <<<<< added
  use mod_sph, only : lmax1  ! <<<<< added

  implicit none
  real(c_double), intent(in) :: dtime
  complex(c_double_complex), intent(in) ::   orb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: h0orb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: h1orb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) ::  gorb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: v2orb(1:(nrad-1), 0:lmax1, 1:nfun)

  integer(c_int), external :: util_omp_nproc
  integer(c_int), external :: util_omp_iproc
  complex(c_double_complex) :: tmp, h012, tfac
  complex(c_double_complex), allocatable :: xmatp(:,:,:)
  integer(c_int) :: nproc, iproc, ifun, jfun, l, llr, ulr, irad, jrad, isub, iact, jact ! <<<<<<< changed


  if (nact == 0) return

  nproc = util_omp_nproc()
  allocate(xmatp(1:nfun, 1:nfun, 0:(nproc-1)))

  !$omp parallel default(shared) private(iproc, llr, ulr, tmp, h012, ifun, jfun, jrad)
  !###########################
  iproc = util_omp_iproc()
  xmatp(1:nfun, 1:nfun, iproc) = czero
  call util_omp_disp(1, irad_ecs-1, llr, ulr)

  do isub = 1, nsub
     do iact = lorb_sub(1, isub), lorb_sub(2, isub)
        ifun = ncore + iact
        do jact = lorb_sub(1, isub), lorb_sub(2, isub)
           jfun = ncore + jact
           if (mval(ifun) == mval(jfun)) then
              tmp = czero
              if (ioorot == 0) then
                 ! <j|h0+h1+h2|i>
                 do l = abs(mval(ifun)), lmax1
                    do irad = llr, ulr
                       h012 = v2orb(irad, l, ifun) + gorb(irad, l, ifun) + h1orb(irad, l, ifun) + h0orb(irad, l, ifun)
                       tmp = tmp + conjg(orb(irad, l, jfun)) * h012
! Sato_ECS
                       h012 = h1orb_out(irad,l,jfun) + h0orb_out(irad,l,jfun)
                       tmp = tmp + conjg(h012) * orb(irad,l,ifun)
! Sato_ECS
                    end do
                 end do
              else if (ioorot == 2) then
                 ! <j|h1+h2|i>
                 do l = abs(mval(ifun)), lmax1
                    do irad = llr, ulr
                       h012 = v2orb(irad, l, ifun) + gorb(irad, l, ifun) + h1orb(irad, l, ifun)
                       tmp = tmp + conjg(orb(irad, l, jfun)) * h012
                    end do
! Sato_ECS
                       h012 = h1orb_out(irad,l,jfun)
                       tmp = tmp + conjg(h012) * orb(irad,l,ifun)
! Sato_ECS
                 end do
              else if (ioorot == 1) then
                 ! <j|h2|i>
                 do l = abs(mval(ifun)), lmax1
                    do irad = llr, ulr
                       h012 = v2orb(irad, l, ifun) + gorb(irad, l, ifun)
                       tmp = tmp + conjg(orb(irad, l, jfun)) * h012
                    end do
                 end do
              end if
              xmatp(jfun, ifun, iproc) = xmatp(jfun, ifun, iproc) + tmp
           end if
        end do
     end do
  end do
  !###########################
  !$omp end parallel

  do iproc = 0, nproc - 1
     do isub = 1, nsub
        do iact = lorb_sub(1, isub), lorb_sub(2, isub)
           ifun = ncore + iact
           do jact = lorb_sub(1, isub), lorb_sub(2, isub)
              jfun = ncore + jact
              if (mval(ifun) == mval(jfun)) then
                 xmat(jfun, ifun) = xmat(jfun, ifun) + xmatp(jfun, ifun, iproc)
              end if
           end do
        end do
     end do
  end do

  do isub = 1, nsub
     do iact = lorb_sub(1, isub), lorb_sub(2, isub)
        ifun = ncore + iact
        do jact = lorb_sub(1, isub), lorb_sub(2, isub)
           jfun = ncore + jact
           if (mval(ifun) == mval(jfun)) then
              fmat(jfun, ifun) = xmat(jfun, ifun)
           end if
        end do
     end do
  end do

  if (icomp == 0) then
     tfac = runit * dtime
  else
     tfac = iunit * dtime
  end if
  do isub = 1, nsub
     do iact = lorb_sub(1, isub), lorb_sub(2, isub)
        ifun = ncore + iact
        do jact = lorb_sub(1, isub), lorb_sub(2, isub)
           jfun = ncore + jact
           if (mval(ifun) == mval(jfun)) then
              xmat(jfun, ifun) = xmat(jfun, ifun) * tfac
           end if
        end do
     end do
  end do

  deallocate(xmatp)

end subroutine hprod_htot_mkxmat_aa1_ecs
!######################################################################
subroutine hprod_htot_mkxmat_aa2_ecs(dtime, orb, h0orb, h1orb, gorb, v2orb, cic, dcic)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : nbas, mval
  use mod_const, only : czero, runit, iunit
  use mod_ormas, only : ncore, nact, nfun, nsub, lorb_sub
  use mod_control, only : icomp, ioorot
  use mod_hprod, only : xmat, fmat, bmat, den1, den2, int1e, int2e, ene_act, h0orb_out, h1orb_out
  use mod_rad, only : ecs_flag, irad_ecs, nrad ! <<<<< added
  use mod_sph, only : lmax1  ! <<<<< added

  implicit none
  real(c_double), intent(in) :: dtime
  complex(c_double_complex), intent(in) ::   orb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: h0orb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: h1orb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) ::  gorb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: v2orb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: cic(1:*)
  complex(c_double_complex), intent(inout) :: dcic(1:*)

  integer(c_int), external :: util_omp_nproc
  integer(c_int), external :: util_omp_iproc
  complex(c_double_complex) :: tmp1, tmp2, tfac
  complex(c_double_complex), allocatable :: xmatp(:,:,:)
  integer(c_int) :: nproc, iproc, ifun, jfun, kfun, l, irad, llr, ulr, isub, jsub, iact, jact, kact

  if (nact == 0) return

  nproc = util_omp_nproc()
  allocate(xmatp(1:nfun, 1:nfun, 0:(nproc-1)))

  !$omp parallel default(shared) private(iproc, llr, ulr, tmp1, tmp2, ifun, jfun)
  !###########################
  iproc = util_omp_iproc()
  xmatp(1:nfun, 1:nfun, iproc) = czero
  call util_omp_disp(1, irad_ecs, llr, ulr)

  do isub = 1, nsub
     do jsub = isub + 1, nsub
        do iact = lorb_sub(1, isub), lorb_sub(2, isub)
           ifun = ncore + iact
           do jact = lorb_sub(1, jsub), lorb_sub(2, jsub)
              jfun = ncore + jact
              if (mval(ifun) == mval(jfun)) then
                 tmp1 = czero
                 tmp2 = czero
                 if (icomp == 1) then
                    do l = abs(mval(ifun)), lmax1
                       do irad = llr, ulr
                          ! Hermitian part is included
                          tmp1 = tmp1 + conjg(orb(irad, l, jfun)) &
                               * (v2orb(irad, l, ifun) + gorb(irad, l, ifun) + h1orb(irad, l, ifun) + h0orb(irad, l, ifun))
                          tmp2 = tmp2 + conjg(orb(irad, l, ifun)) &
                               * (v2orb(irad, l, jfun) + gorb(irad, l, jfun) + h1orb(irad, l, jfun) + h0orb(irad, l, jfun))
! Sato_ECS
                          tmp1 = tmp1 + conjg(h1orb_out(irad,l,jfun)+h0orb_out(irad,l,jfun)) * orb(irad,l,ifun)
! Bug fixed:2018/02/07    tmp1 = tmp1 + conjg(h1orb_out(irad,l,ifun)+h0orb_out(irad,l,ifun)) * orb(irad,l,jfun)
                          tmp2 = tmp2 + conjg(h1orb_out(irad,l,ifun)+h0orb_out(irad,l,ifun)) * orb(irad,l,jfun)
! Sato_ECS
                       end do
                    end do
                 else
                    ! Hermitian part is included
                    do l = abs(mval(ifun)), lmax1
                       do irad = llr, ulr
                          tmp1 = tmp1 + conjg(orb(irad, l, jfun)) &
                               * (v2orb(irad, l, ifun) + gorb(irad, l, ifun) + h1orb(irad, l, ifun) + h0orb(irad, l, ifun))
                          tmp2 = tmp2 + conjg(orb(irad, l, ifun)) &
                               * (v2orb(irad, l, jfun) + gorb(irad, l, jfun) + h1orb(irad, l, jfun) + h0orb(irad, l, jfun))
                       end do
                    end do
                 end if
                 xmatp(jfun, ifun, iproc) = xmatp(jfun, ifun, iproc) + tmp1
                 xmatp(ifun, jfun, iproc) = xmatp(ifun, jfun, iproc) + tmp2
              end if
           end do
        end do
     end do
  end do
  !###########################
  !$omp end parallel

  do iproc = 0, nproc - 1
     do isub = 1, nsub
        do jsub = isub + 1, nsub
           do iact = lorb_sub(1, isub), lorb_sub(2, isub)
              ifun = ncore + iact
              do jact = lorb_sub(1, jsub), lorb_sub(2, jsub)
                 jfun = ncore + jact
                 if (mval(ifun) == mval(jfun)) then
                    xmat(jfun, ifun) = xmat(jfun, ifun) + xmatp(jfun, ifun, iproc)
                    xmat(ifun, jfun) = xmat(ifun, jfun) + xmatp(ifun, jfun, iproc)
                 end if
              end do
           end do
        end do
     end do
  end do

  do isub = 1, nsub
     do jsub = isub + 1, nsub
        do iact = lorb_sub(1, isub), lorb_sub(2, isub)
           ifun = ncore + iact
           do jact = lorb_sub(1, jsub), lorb_sub(2, jsub)
              jfun = ncore + jact
              if (mval(ifun) == mval(jfun)) then
                 fmat(jfun, ifun) = xmat(jfun, ifun)
                 fmat(ifun, jfun) = xmat(ifun, jfun)
              end if
           end do
        end do
     end do
  end do

  ! TD-GBT matrix
  do isub = 1, nsub
     do jsub = isub + 1, nsub
        do iact = lorb_sub(1, isub), lorb_sub(2, isub)
           ifun = ncore + iact
           do jact = lorb_sub(1, jsub), lorb_sub(2, jsub)
              jfun = ncore + jact
              if (mval(ifun) == mval(jfun)) then
                 tmp1 = czero
                 tmp2 = czero
                 do kact = 1, nact
                    kfun = ncore + kact
                    if (mval(jfun) == mval(kfun)) then
                       tmp1 = tmp1 + fmat(jfun, kfun) * den1(kact, iact)
                       tmp2 = tmp2 + fmat(ifun, kfun) * den1(kact, jact)
                    end if
                 end do
                 bmat(jfun, ifun) = tmp1 - conjg(tmp2)
                 bmat(ifun, jfun) = -conjg(bmat(jfun, ifun))
              else
                 bmat(jfun, ifun) = czero
                 bmat(ifun, jfun) = czero
              end if
           end do
        end do
     end do
  end do
!DEBUG
!  write(6, "('# hprod_htot_mkxmat_aa2: fmat')")
!  call util_matoutc(nfun, fmat)
!  write(6, "('# hprod_htot_mkxmat_aa2: bmat')")
!  call util_matoutc(nfun, bmat)
!DEBUG

  !old call ormas_xact2(dtime, cic, den1, den2, bmat, dcic)
  !old call ormas_xact2ras(dtime, cic, den1, den2, bmat, dcic)
  call ormas_xact2(dtime, int1e, int2e, ene_act, cic, den1, den2, bmat, dcic)

  ! P + NR contribution
  if (icomp == 0) then
     tfac = runit * dtime
  else
     tfac = iunit * dtime
  end if

  do isub = 1, nsub
     do jsub = isub + 1, nsub
        do iact = lorb_sub(1, isub), lorb_sub(2, isub)
           ifun = ncore + iact
           do jact = lorb_sub(1, jsub), lorb_sub(2, jsub)
              jfun = ncore + jact
              if (mval(ifun) == mval(jfun)) then
                 xmat(jfun, ifun) = xmat(jfun, ifun) * tfac + bmat(jfun, ifun)
                 xmat(ifun, jfun) = xmat(ifun, jfun) * tfac + bmat(ifun, jfun)
              end if
           end do
        end do
     end do
  end do

  deallocate(xmatp)

end subroutine hprod_htot_mkxmat_aa2_ecs
!######################################################################
subroutine hprod_htot_mkxmat_ca_ecs(dtime, orb, h0orb, h1orb, gorb, v2orb)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : nbas, mval
  use mod_control, only : icomp, ioorot
  use mod_const, only : czero, runit, iunit
  use mod_ormas, only : nfcore, ndcore, ncore, nact, nfun, nelact
  use mod_hprod, only : xmat, fmat, bmat, den1, rrden, h0orb_out, h1orb_out

  use mod_rad, only : irad_ecs, nrad ! <<<<< added
  use mod_sph, only : lmax1  ! <<<<< added

  implicit none
  real(c_double), intent(in) :: dtime
  complex(c_double_complex), intent(in) ::   orb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: h0orb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: h1orb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) ::  gorb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: v2orb(1:(nrad-1), 0:lmax1, 1:nfun)

  integer(c_int), external :: util_omp_nproc
  integer(c_int), external :: util_omp_iproc
  complex(c_double_complex) :: h012C, h012A, tmpAC, tmpCA, tfac
  complex(c_double_complex), allocatable :: fmatp(:,:,:)
  integer(c_int) :: nproc, iproc, ifun, jfun, kfun, jact, kact, l, irad, llr, ulr

  ! trivial case
  if (nact * 2 == nelact(3)) then
     call hprod_htot_mkxmat_ca0_ecs(dtime, orb, h0orb, h1orb, gorb, v2orb)
     return
  end if

  if (ncore == 0 .or. nact == 0) return
  fmat((nfcore+1):ncore, (ncore+1):nfun) = czero
  fmat((ncore+1):nfun, (nfcore+1):ncore) = czero
  bmat((nfcore+1):ncore, (ncore+1):nfun) = czero
  bmat((ncore+1):nfun, (nfcore+1):ncore) = czero
  xmat((nfcore+1):ncore, (ncore+1):nfun) = czero
  xmat((ncore+1):nfun, (nfcore+1):ncore) = czero

  if (icomp == 0) then
     tfac = runit * dtime
  else
     tfac = iunit * dtime
  end if

  nproc = util_omp_nproc()
  allocate(fmatp(1:nfun, 1:nfun, 0:(nproc-1)))

!  if (ndcore > 0 .and. icomp == 1) then
!     write(6,"('hprod_htot_mkxmat_ca: hermitian part is included.')")
!!     write(6,"('hprod_htot_mkxmat_ca: hermitian part is excluded.')")
!  end if

  !$omp parallel default(shared) private(iproc, llr, ulr, h012C, h012A, tmpAC, tmpCA)
  !###########################
  iproc = util_omp_iproc()
  call util_omp_disp(1, irad_ecs-1, llr, ulr)

  do ifun = nfcore + 1, ncore
     do jfun = ncore + 1, nfun
        tmpAC = czero
        tmpCA = czero
        if (mval(jfun) == mval(ifun)) then
           ! <t|F_i|i>
           ! <i|F_t|t>
           if (icomp == 1) then
              do l = abs(mval(ifun)), lmax1
                 do irad = llr, ulr                    
!WATCH HERE! >>>>
                    ! Hermitian part is canceled out.
                    ! h012C = v2orb(ibas, ifun)
                    ! h012A = v2orb(ibas, jfun)
                    ! Hermitian part is included
                 h012C = v2orb(irad, l, ifun) + gorb(irad, l, ifun) + h1orb(irad, l, ifun) + h0orb(irad, l, ifun)
                 h012A = v2orb(irad, l, jfun) + gorb(irad, l, jfun) + h1orb(irad, l, jfun) + h0orb(irad, l, jfun)
!WATCH HERE! <<<<
                 tmpAC = tmpAC + conjg(orb(irad, l, jfun)) * h012C
                 tmpCA = tmpCA + conjg(orb(irad, l, ifun)) * h012A
! Sato_ECS
                 h012C = h1orb_out(irad,l,jfun) + h0orb_out(irad,l,jfun)
                 h012A = h1orb_out(irad,l,ifun) + h0orb_out(irad,l,ifun)
                 tmpAC = tmpAC + conjg(h012C) * orb(irad,l,ifun)
                 tmpCA = tmpCA + conjg(h012A) * orb(irad,l,jfun)
! Sato_ECS
              end do
           end do
           else
              ! Hermitian part is included
              do l = abs(mval(ifun)), lmax1
                 do irad = llr, ulr
                    h012C = v2orb(irad, l, ifun) + gorb(irad, l, ifun) + h1orb(irad, l, ifun) + h0orb(irad, l, ifun)
                    h012A = v2orb(irad, l, jfun) + gorb(irad, l, jfun) + h1orb(irad, l, jfun) + h0orb(irad, l, jfun)
                    tmpAC = tmpAC + conjg(orb(irad, l, jfun)) * h012C
                    tmpCA = tmpCA + conjg(orb(irad, l, ifun)) * h012A
              end do
           end do
           end if
        end if
        fmatp(jfun, ifun, iproc) = tmpAC
        fmatp(ifun, jfun, iproc) = tmpCA
     end do
  end do
  !###########################
  !$omp end parallel

  do iproc = 0, nproc - 1
     do ifun = nfcore + 1, ncore
        do jfun = ncore + 1, nfun
           fmat(jfun, ifun) = fmat(jfun, ifun) + fmatp(jfun, ifun, iproc)
           fmat(ifun, jfun) = fmat(ifun, jfun) + fmatp(ifun, jfun, iproc)
        end do
     end do
  end do

  ! P contribution
  do ifun = nfcore + 1, ncore
     do jfun = ncore + 1, nfun
        xmat(jfun, ifun) = fmat(jfun, ifun) * tfac
        xmat(ifun, jfun) = fmat(ifun, jfun) * tfac
     end do
  end do

  ! F_AC = 2F_AC - Den1 F_CA^*
  do ifun = nfcore + 1, ncore
     do jact = 1, nact
        jfun = ncore + jact
        fmat(jfun, ifun) = fmat(jfun, ifun) + fmat(jfun, ifun)
        tmpCA = czero
        do kact = 1, nact
           kfun = ncore + kact
           tmpCA = tmpCA + den1(jact, kact) * conjg(fmat(ifun, kfun))
        end do
        fmat(jfun, ifun) = fmat(jfun, ifun) - tmpCA
     end do
  end do

  ! B = 1/(2 - Den1)(2F_AC - Den1 F_CA^*)
  do ifun = nfcore + 1, ncore
     do jact = 1, nact
        jfun = ncore + jact
        tmpCA = czero
        do kact = 1, nact
           kfun = ncore + kact
           tmpCA = tmpCA + rrden(jact, kact) * fmat(kfun, ifun)
        end do
        bmat(jfun, ifun) = - tfac * tmpCA
        bmat(ifun, jfun) = - conjg(bmat(jfun, ifun))
     end do
  end do

  ! R contribution
  do ifun = nfcore + 1, ncore
     do jfun = ncore + 1, nfun
        xmat(jfun, ifun) = xmat(jfun, ifun) + bmat(jfun, ifun)
        xmat(ifun, jfun) = xmat(ifun, jfun) + bmat(ifun, jfun)
!        xmat(jfun, ifun) = xmat(jfun, ifun) + bmat(jfun, ifun)
!        xmat(ifun, jfun) = xmat(ifun, jfun) - bmat(ifun, jfun)
!        xmat(jfun, ifun) = xmat(jfun, ifun) - bmat(jfun, ifun)
!        xmat(ifun, jfun) = xmat(ifun, jfun) - bmat(ifun, jfun)
     end do
  end do

!DEBUG
!  write(6, "('x-CA:')")
!  do ifun = nfcore + 1, ncore
!     do jact = 1, nact
!        jfun = ncore + jact
!        write(6, "(2f20.10)", advance = 'no') dble(fmat(jfun, ifun))
!     end do
!     do jact = 1, nact
!        jfun = ncore + jact
!        write(6, "(2f20.10)", advance = 'no') aimag(fmat(jfun, ifun))
!     end do
!     do jact = 1, nact
!        jfun = ncore + jact
!        write(6, "(2f20.10)", advance = 'no') dble(bmat(jfun, ifun))
!     end do
!     do jact = 1, nact
!        jfun = ncore + jact
!        write(6, "(2f20.10)", advance = 'no') aimag(bmat(jfun, ifun))
!     end do
!     do jact = 1, nact
!        jfun = ncore + jact
!        write(6, "(2f20.10)", advance = 'no') dble(xmat(jfun, ifun))
!     end do
!     do jact = 1, nact
!        jfun = ncore + jact
!        write(6, "(2f20.10)", advance = 'no') aimag(xmat(jfun, ifun))
!     end do
!     write(6, *)
!  end do
!DEBUG

  deallocate(fmatp)

end subroutine hprod_htot_mkxmat_ca_ecs
!######################################################################
subroutine hprod_htot_mkxmat_ca0_ecs(dtime, orb, h0orb, h1orb, gorb, v2orb)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : nbas, mval
  use mod_const, only : czero, runit, iunit
  use mod_ormas, only : nfcore, ncore, nact, nfun
  use mod_control, only : icomp, ioorot
  use mod_hprod, only : xmat, h0orb_out, h1orb_out

  use mod_rad, only : irad_ecs, nrad ! <<<<< added
  use mod_sph, only : lmax1  ! <<<<< added

  implicit none
  real(c_double), intent(in) :: dtime
  complex(c_double_complex), intent(in) ::   orb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: h0orb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: h1orb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) ::  gorb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: v2orb(1:(nrad-1), 0:lmax1, 1:nfun)

  integer(c_int), external :: util_omp_nproc
  integer(c_int), external :: util_omp_iproc
  complex(c_double_complex) :: h012C, h012A, tmpAC, tmpCA, tfac
  complex(c_double_complex), allocatable :: xmatp(:,:,:)
  integer(c_int) :: nproc, iproc, ifun, jfun, kfun, jact, kact, l, irad, llr, ulr

  if (ncore == nfcore) return
  if (ncore == 0) return
  if (nact == 0) return

  xmat((nfcore+1):ncore, (ncore+1):nfun) = czero
  xmat((ncore+1):nfun, (nfcore+1):ncore) = czero

  if (icomp == 0) then
     tfac = runit * dtime
  else
     tfac = iunit * dtime
  end if

  nproc = util_omp_nproc()
  allocate(xmatp(1:nfun, 1:nfun, 0:(nproc-1)))

  !$omp parallel default(shared) private(iproc, llr, ulr, h012C, h012A, tmpAC, tmpCA)
  !###########################
  iproc = util_omp_iproc()
  call util_omp_disp(1, irad_ecs-1, llr, ulr)

  do ifun = nfcore + 1, ncore
     do jfun = ncore + 1, nfun
        tmpAC = czero
        tmpCA = czero
        if (mval(jfun) == mval(ifun)) then
           if (ioorot == 0) then
              ! <j|h0+h1+h2|i>
              do l = abs(mval(ifun)), lmax1
                 do irad = llr, ulr
                    h012C = v2orb(irad, l, ifun) + gorb(irad, l, ifun) + h1orb(irad, l, ifun) + h0orb(irad, l, ifun)
                    h012A = v2orb(irad, l, jfun) + gorb(irad, l, jfun) + h1orb(irad, l, jfun) + h0orb(irad, l, jfun)
                    tmpAC = tmpAC + conjg(orb(irad, l, jfun)) * h012C
                    tmpCA = tmpCA + conjg(orb(irad, l, ifun)) * h012A
! Sato_ECS
                    h012C = h1orb_out(irad,l,jfun) + h0orb_out(irad,l,jfun)
                    h012A = h1orb_out(irad,l,ifun) + h0orb_out(irad,l,ifun)
                    tmpAC = tmpAC + conjg(h012C) * orb(irad,l,ifun)
                    tmpCA = tmpCA + conjg(h012A) * orb(irad,l,jfun)
! Sato_ECS
                 end do
              end do
           else if (ioorot == 2) then
              ! <j|h1+h2|i>
              do l = abs(mval(ifun)), lmax1
                 do irad = llr, ulr
                    h012C = v2orb(irad, l, ifun) + gorb(irad, l, ifun) + h1orb(irad, l, ifun)
                    h012A = v2orb(irad, l, jfun) + gorb(irad, l, jfun) + h1orb(irad, l, jfun)
                    tmpAC = tmpAC + conjg(orb(irad, l, jfun)) * h012C
                    tmpCA = tmpCA + conjg(orb(irad, l, ifun)) * h012A
! Sato_ECS
                    h012C = h1orb_out(irad,l,jfun)
                    h012A = h1orb_out(irad,l,ifun)
                    tmpAC = tmpAC + conjg(h012C) * orb(irad,l,ifun)
                    tmpCA = tmpCA + conjg(h012A) * orb(irad,l,jfun)
! Sato_ECS
              end do
           end do
           else if (ioorot == 1) then
              ! <j|h2|i>
              do l = abs(mval(ifun)), lmax1
                 do irad = llr, ulr
                    h012C = v2orb(irad, l, ifun) + gorb(irad, l, ifun)
                    h012A = v2orb(irad, l, jfun) + gorb(irad, l, jfun)
                    tmpAC = tmpAC + conjg(orb(irad, l, jfun)) * h012C
                    tmpCA = tmpCA + conjg(orb(irad, l, ifun)) * h012A
              end do
           end do
           end if
        end if
        xmatp(jfun, ifun, iproc) = tmpAC
        xmatp(ifun, jfun, iproc) = tmpCA
     end do
  end do
  !###########################
  !$omp end parallel

  do iproc = 0, nproc - 1
     do ifun = nfcore + 1, ncore
        do jfun = ncore + 1, nfun
           xmat (jfun, ifun) = xmat (jfun, ifun) + xmatp(jfun, ifun, iproc)
           xmat (ifun, jfun) = xmat (ifun, jfun) + xmatp(ifun, jfun, iproc)
        end do
     end do
  end do

  do ifun = nfcore + 1, ncore
     do jfun = ncore + 1, nfun
        xmat(jfun, ifun) = xmat(jfun, ifun) * tfac
        xmat(ifun, jfun) = xmat(ifun, jfun) * tfac
     end do
  end do

  deallocate(xmatp)

end subroutine hprod_htot_mkxmat_ca0_ecs
!######################################################################

