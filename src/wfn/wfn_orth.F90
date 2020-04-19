!///////////////////////////////////////////////////////////////////////
subroutine wfn_orth_cic(cic)

  use, intrinsic :: iso_c_binding
  use mod_const, only : zero, one
  use mod_ormas, only : lcic, tdcc

  implicit none
  complex(c_double_complex), intent(inout) :: cic(1:lcic)

  real(c_double) :: norm
  integer(c_long) :: idet

  if (tdcc) then
!     write(6,"('wfn_orth_cic: cc amplitudes left unnormalized.')")
     return
  end if
!     write(6,"('WARNING: wfn_orth_cic: cc amplitude normalized.')")

  norm = zero
  do idet = 1, lcic
     norm = norm + dble(conjg(cic(idet)) * cic(idet))
  end do
  norm = one / sqrt(norm)
  do idet = 1, lcic
     cic(idet) = cic(idet) * norm
  end do

end subroutine wfn_orth_cic
!///////////////////////////////////////////////////////////////////////
subroutine wfn_orth_orb(orb)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad, wrad, ecs_flag, irad_ecs, cwrad
  use mod_sph, only : lmax1
  use mod_bas, only : mval
  use mod_ormas, only : nfcore, nfun
  use mod_const, only : zero, one, czero
  use mod_control, only : fedvr_normalized

  implicit none
  complex(c_double_complex), intent(inout) :: orb(1:(nrad-1), 0:lmax1, 1:nfun)

  real(c_double) :: norm
  complex(c_double_complex) :: tmp, ovlp
  integer(c_long) :: ifun, jfun, l, mi, mj, irad
  real(c_double), allocatable :: bra_wrad(:)

!DEBUG
! symmetric orthonormalization X = S^{-0.5}.
! call wfn_orth_orb_symm(nfun, nbas, orb)
! return
!DEBUG

  allocate(bra_wrad(1:nrad-1))

  if(ecs_flag == 0) then
     
!    do ifun = 1, nfun
     do ifun = nfcore + 1, nfun
        mi = mval(ifun)

        do jfun = 1, ifun - 1
           mj = mval(jfun)
           if (mj == mi) then
              
              tmp = czero
              if (fedvr_normalized) then
                 do l = abs(mi), lmax1
                    do irad = 1, nrad - 1
                       tmp = tmp + conjg(orb(irad, l, jfun)) * orb(irad, l, ifun)
                    end do
                 end do
              else
                 do l = abs(mi), lmax1
                    do irad = 1, nrad - 1
                       tmp = tmp + conjg(orb(irad, l, jfun)) * orb(irad, l, ifun) * wrad(irad)
                    end do
                 end do
              end if
              ovlp = tmp
              !debug
              !        write(6, "('wfn_orth_orb: ovlp: ', 2i5, 2f20.10)") jfun, ifun, ovlp
              !debug
              do l = abs(mi), lmax1
                 do irad = 1, nrad - 1
                    orb(irad, l, ifun) = orb(irad, l, ifun) - orb(irad, l, jfun) * ovlp
                 end do
              end do
           end if
        end do
        
        norm = zero
        if (fedvr_normalized) then
           do l = abs(mi), lmax1
              do irad = 1, nrad - 1
                 norm = norm + dble(conjg(orb(irad, l, ifun)) * orb(irad, l, ifun))
              end do
           end do
        else
           do l = abs(mi), lmax1
              do irad = 1, nrad - 1
                 norm = norm + dble(conjg(orb(irad, l, ifun)) * orb(irad, l, ifun)) * wrad(irad)
              end do
           end do
        end if
        !debug:
        !     write(6, "('norm: ', 2i5, 2f20.10)") ifun, ifun, sqrt(norm)
        !debug
        norm = one / sqrt(norm)
        do l = abs(mi), lmax1
           do irad = 1, nrad - 1
              orb(irad, l, ifun) = orb(irad, l, ifun) * norm
           end do
        end do
     end do

  else if (ecs_flag == 1) then

     do irad = 1, nrad - 1
        bra_wrad(irad) = wrad(irad) / dble(sqrt(conjg(cwrad(irad)) * cwrad(irad))) 
        !debug
        !write(*,'(f20.10)') bra_wrad(irad)
        !debug
     end do
     
     do ifun = 1, nfun
        mi = mval(ifun)

        do jfun = 1, ifun - 1
           mj = mval(jfun)
           if (mj == mi) then

              tmp = czero
              do l = abs(mi), lmax1
                 do irad = 1, nrad - 1
                    tmp = tmp + conjg(orb(irad, l, jfun)) * bra_wrad(irad) * orb(irad, l, ifun)
                 end do
                 

              end do
              ovlp = tmp
              !debug        write(6, "('wfn_orth_orb: ovlp: ', 2i5, 2f20.10)") jfun, ifun, ovlp
              do l = abs(mi), lmax1
                 do irad = 1, nrad - 1
                    orb(irad, l, ifun) = orb(irad, l, ifun) - orb(irad, l, jfun) * ovlp
                 end do
              end do
           end if
        end do

        norm = zero
        do l = abs(mi), lmax1
           do irad = 1, nrad - 1
              norm = norm + dble(conjg(orb(irad, l, ifun)) * bra_wrad(irad) * orb(irad, l, ifun))
           end do
        end do
        !debug:
        !     write(6, "('norm: ', 2i5, 2f20.10)") ifun, ifun, sqrt(norm)
        !debug
        norm = one / sqrt(norm)
        do l = abs(mi), lmax1
           do irad = 1, nrad - 1
              orb(irad, l, ifun) = orb(irad, l, ifun) * norm
           end do
        end do
     end do

!DEBUG
     write(*,'(f20.10)') norm
!DEBUG     
  end if

  deallocate(bra_wrad)

end subroutine wfn_orth_orb
!///////////////////////////////////////////////////////////////////////
subroutine wfn_orth_orb_symm(orb)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : nbas
  use mod_ormas, only : nfun

  implicit none
  complex(c_double_complex), intent(inout) :: orb(1:nbas, 1:nfun)

  complex(c_double_complex) :: tmp
  integer(c_long) :: ifun, jfun, kfun, ibas

  complex(c_double_complex), allocatable :: seig(:)
  complex(c_double_complex), allocatable :: smat(:,:)
  complex(c_double_complex), allocatable :: umat(:,:)
  complex(c_double_complex), allocatable :: torb(:,:)

  write(6, "('wfn_orth_orb_symm may include bug!')")
  stop

  allocate(seig(1:nfun))
  allocate(smat(1:nfun, 1:nfun))
  allocate(umat(1:nfun, 1:nfun))
  allocate(torb(1:nbas, 1:nfun))
  seig(1:nfun) = (0.d+0, 0.d+0)
  smat(1:nfun, 1:nfun) = (0.d+0, 0.d+0)
  umat(1:nfun, 1:nfun) = (0.d+0, 0.d+0)

  ! DEBUG
  ! write(6, "('before symmetric orthonormalization:')")
  ! DEBUG
  call wfn_orth_smat(orb, smat)

!  do ifun = 1, nfun
!     do jfun = 1, nfun
!        tmp = (0.d+0, 0.d+0)
!        do ibas = 1, nbas
!           tmp = tmp + conjg(orb(ibas, jfun)) * orb(ibas, ifun)
!        end do
!        smat(jfun, ifun) = tmp
!     end do
!  end do
!
!  !DEBUG
!  write(6, "('smat-real-1')")
!  do ifun = 1, nfun
!     do jfun = 1, nfun
!        write(6, "(f10.5)", advance = 'no') dble(smat(jfun, ifun))
!     end do
!     write(6, *)
!  end do
!  write(6, "('smat-imag-1')")
!  do ifun = 1, nfun
!     do jfun = 1, nfun
!        write(6, "(f10.5)", advance = 'no') aimag(smat(jfun, ifun))
!     end do
!     write(6, *)
!  end do
!  !stop 'stop @ wfn_orth for debug...'
!  !DEBUG

  call lapack_zheev_dsc(nfun, smat, umat)
  ! DEBUG
  ! write(6, "('eigen-solutions of smat (R):')")
  ! do kfun = 1, nfun
  !    write(6, "(5x, i12)", advance = 'no') kfun
  ! end do
  ! write(6, *)
  ! do kfun = 1, nfun
  !    write(6, "(5x, f12.5)", advance = 'no') dble(smat(kfun, kfun))
  ! end do
  ! write(6, *)
  ! do ifun = 1, nfun
  !    do kfun = 1, nfun
  !       write(6, "(i5, f12.5)", advance = 'no') ifun, dble(umat(ifun, kfun))
  !    end do
  !    write(6, *)
  ! end do
  ! write(6, "('eigen-solutions of smat (I):')")
  ! do kfun = 1, nfun
  !    write(6, "(5x, i12)", advance = 'no') kfun
  ! end do
  ! write(6, *)
  ! do kfun = 1, nfun
  !    write(6, "(5x, f12.5)", advance = 'no') aimag(smat(kfun, kfun))
  ! end do
  ! write(6, *)
  ! do ifun = 1, nfun
  !    do kfun = 1, nfun
  !       write(6, "(i5, f12.5)", advance = 'no') ifun, aimag(umat(ifun, kfun))
  !    end do
  !    write(6, *)
  ! end do
  ! DEBUG

  do kfun = 1, nfun
     seig(kfun) = smat(kfun, kfun)
  end do

  smat(1:nfun, 1:nfun) = (0.d+0, 0.d+0)
  do kfun = 1, nfun
     tmp = seig(kfun)
     if (abs(tmp) > 1.D-10) then
        tmp = tmp ** (-0.5D+0)
        do ifun = 1, nfun
           do jfun = 1, nfun
              smat(jfun, ifun) = smat(jfun, ifun) &
                             & + umat(jfun, kfun) * tmp &
                       & * conjg(umat(ifun, kfun))
           end do
        end do
     end if
  end do

  do ifun = 1, nfun
     do ibas = 1, nbas
        torb(ibas, ifun) = orb(ibas, ifun)
        orb (ibas, ifun) = (0.d+0, 0.d+0)
     end do
  end do

  do ifun = 1, nfun
     do jfun = 1, nfun
        do ibas = 1, nbas
           orb(ibas, ifun) = orb(ibas, ifun) + torb(ibas, jfun) * smat(jfun, ifun)
        end do
     end do
  end do

  ! DEBUG
  ! write(6, "('after symmetric orthonormalization:')")
  ! call wfn_orth_smat(orb, smat)
  ! DEBUG

  deallocate(umat)
  deallocate(smat)
  deallocate(seig)
  deallocate(torb)

end subroutine wfn_orth_orb_symm
!///////////////////////////////////////////////////////////////////////
subroutine wfn_orth_smat(orb, smat)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : nbas
  use mod_ormas, only : nfun

  implicit none
  complex(c_double_complex), intent(in) :: orb(1:nbas, 1:nfun)
  complex(c_double_complex), intent(out) :: smat(1:nfun, 1:nfun)
  
  complex(c_double_complex) :: tmp
  integer(c_long) :: ifun, jfun, ibas

  write(6, "('wfn_orth_smat may include bug!')")
  stop

  smat(1:nfun, 1:nfun) = (0.d+0, 0.d+0)

  do ifun = 1, nfun
     do jfun = 1, nfun
        tmp = (0.d+0, 0.d+0)
        do ibas = 1, nbas
           tmp = tmp + conjg(orb(ibas, jfun)) * orb(ibas, ifun)
        end do
        smat(jfun, ifun) = tmp
     end do
  end do

  !DEBUG
  write(6, "('smat-real')")
  do ifun = 1, nfun
     do jfun = 1, nfun
        write(6, "(f10.5)", advance = 'no') dble(smat(jfun, ifun))
     end do
     write(6, *)
  end do
  write(6, "('smat-imag')")
  do ifun = 1, nfun
     do jfun = 1, nfun
        write(6, "(f10.5)", advance = 'no') aimag(smat(jfun, ifun))
     end do
     write(6, *)
  end do
  !stop 'stop @ wfn_orth for debug...'
  !DEBUG

end subroutine wfn_orth_smat
!///////////////////////////////////////////////////////////////////////
