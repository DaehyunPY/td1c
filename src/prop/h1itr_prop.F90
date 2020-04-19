!######################################################################
subroutine h1itr_prop(icomp, isplit, igauge, iprojfc, lfield, dtime, maxcyc, &
                      thresh, cnpiv, cninv, wfn)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : nbas
  use mod_sph, only : lmax1
  use mod_rad, only : nrad, wrad
  use mod_ormas, only : nfcore, nfun
  use mod_const, only : zero, half, czero, runit, iunit
  use mod_control, only : fedvr_normalized

  implicit none
  integer(c_long), intent(in) :: icomp, isplit, igauge, iprojfc, maxcyc, cnpiv(1:*)
  real(c_double), intent(in) :: lfield(9), dtime, thresh
  complex(c_double_complex), intent(in) :: cninv(1:*)
  complex(c_double_complex), intent(inout) :: wfn(1:(nrad-1), 0:lmax1, 1:nfun)

  complex(c_double_complex) :: zfield, zfac1, zfac2, tmp
  integer(c_long) :: ncyc, icyc, ifun, irad, l, size1, sizef, sized
  logical, parameter :: debug = .true.
  real(c_double) :: rmsres, maxres
  logical projfc
  complex(c_double_complex), allocatable :: orb(:,:,:)
  complex(c_double_complex), allocatable :: horb(:,:,:)
  complex(c_double_complex), allocatable :: torb(:,:,:)

  !DEBUG
  !call h1itr_prop2(icomp, isplit, igauge, iprojfc, lfield, dtime, maxcyc, &
  !                 thresh, cnpiv, cninv, orb, horb, torb, wfn)
  !return
  !DEBUG

  allocate(orb(1:(nrad-1), 0:lmax1, 1:nfun))
  allocate(horb(1:(nrad-1), 0:lmax1, 1:nfun))
  allocate(torb(1:(nrad-1), 0:lmax1, 1:nfun))

  projfc = iprojfc == 1
  size1 = nbas * nfun
  sizef = nbas * nfcore
  sized = size1 - sizef
  zfield = lfield(3)
  if (icomp == 1) then
     zfac1 = - iunit * dtime * half
     zfac2 = - iunit * dtime * half * zfield
  else
     zfac1 = - runit * dtime * half
     zfac2 = - runit * dtime * half * zfield
  end if

  ! initialization
  call zclear_omp(size1, horb)
  call zcopy_omp(size1, wfn, orb)
  call zcopy_omp(size1, wfn, torb)

  !debug if (debug) then
  !debug    write(6, "(i5, ':')", advance = 'no') -2
  !debug    do ifun = 1, nfun
  !debug       call zdotc_omp(nbas, orb(1, ifun), orb(1, ifun), tmp)
  !debug       write(6, "(e20.5)", advance = 'no') sqrt(abs(tmp))
  !debug    end do
  !debug    write(6, *)
  !debug end if

  ! 1 - i * H * dt/2
  call hprod_tprod_dyn(orb, horb)
  if (isplit == 1 .and. igauge == 0) call hprod_zprod_dyn(zfield, orb, horb)
  if (isplit == 1 .and. igauge == 1) call hprod_pzprod_dyn(zfield, orb, horb)

  if (fedvr_normalized) then
     call zaxpy_omp(sized, zfac1, horb(1,0,nfcore+1), orb(1,0,nfcore+1))
  else
     do ifun = nfcore + 1, nfun
        do l = 0, lmax1
           do irad = 1, nrad - 1
              orb(irad, l, ifun) = orb(irad, l, ifun) * wrad(irad) + zfac1 * horb(irad, l, ifun)
           end do
        end do
     end do
  end if

  !debug if (debug) then
  !debug    write(6, "(i5, ':')", advance = 'no') -1
  !debug    do ifun = 1, nfun
  !debug       call zdotc_omp(nbas, orb(1, ifun), orb(1, ifun), tmp)
  !debug       write(6, "(e20.5)", advance = 'no') sqrt(abs(tmp))
  !debug    end do
  !debug    write(6, *)
  !debug end if

  ! 1 / (1 + i * H0 dt /2)
  call h1itr_cninv(cnpiv, cninv, orb)
  call zcopy_omp(size1, orb, torb)

  !debug if (debug) then
  !debug    write(6, "(i5, ':')", advance = 'no') 0
  !debug    do ifun = 1, nfun
  !debug       call zdotc_omp(nbas, orb(1, ifun), orb(1, ifun), tmp)
  !debug       write(6, "(e20.5)", advance = 'no') sqrt(abs(tmp))
  !debug    end do
  !debug    write(6, *)
  !debug end if

  if (isplit == 1) then
     do icyc = 1, maxcyc
        ncyc = icyc
        maxres = zero
        rmsres = zero
        if (debug) write(6, "('h1itr_prop: ', i5, ':')", advance = 'no') icyc
        do ifun = 1, nfun
           if (fedvr_normalized) then
              call zdotc_omp(nbas, orb(1,0,ifun), orb(1,0,ifun), tmp)
           else
              tmp = czero
              do l = 0, lmax1
                 do irad = 1, nrad - 1
                    tmp = tmp + conjg(orb(irad, l, ifun)) * orb(irad, l, ifun) * wrad(irad)
                 end do
              end do
           end if
           if (debug) write(6, "(e20.5)", advance = 'no') sqrt(abs(tmp))

           maxres = max(maxres, sqrt(abs(tmp)))
           rmsres = rmsres + abs(tmp)
        end do
        if (debug) write(6, *)

        rmsres = sqrt(rmsres / nfun)
        if (rmsres < thresh) exit

        ! i * V * dt/2
        call zclear_omp(size1, horb)
        if (projfc) then
           if (igauge == 0) call hprod_zprod_all(zfac2, orb, horb);
           if (igauge == 1) call hprod_pzprod_all(zfac2, orb, horb);
        else
           if (igauge == 0) call hprod_zprod_dyn(zfac2, orb, horb);
           if (igauge == 1) call hprod_pzprod_dyn(zfac2, orb, horb);
        end if

        ! 1 / (1 + i * H0 dt /2)
        call zcopy_omp(size1, horb, orb);
        call h1itr_cninv(cnpiv, cninv, orb);
        call zadd_omp(size1, orb, torb);

     end do
     if (debug) write(6, "('h1itr_prop: norm after ', i5, ' cycles: ' 2e20.5)") icyc, rmsres, maxres
  end if

  ! frozen-core correction
  if (projfc) call h1itr_prop_fcore(wfn, torb)

  call zcopy_omp(sized, torb(1,0,nfcore + 1), wfn(1,0,nfcore + 1))
  !debug if (nfcore > 0 .and. igauge == 1) call h1itr_prop_fcrot(lfield, dtime, wfn)

  !debug if (debug) then
  !debug    write(6, "('orthonormalization check')")
  !debug    do ifun = 1, nfun
  !debug       do jfun = 1, nfun
  !debug          if (mval(jfun) == mval(ifun)) then
  !debug             call zdotc_omp(nbas, wfn(1, ifun), wfn(1, jfun), tmp)
  !debug          else
  !debug             tmp = czero
  !debug          end if
  !debug          write(6, "(2i5, 2e20.5)") ifun, jfun, tmp
  !debug       end do
  !debug    end do
  !debug end if

  deallocate(orb)
  deallocate(horb)
  deallocate(torb)
!  stop "@ h1itr_prop..."

end subroutine h1itr_prop
!######################################################################
subroutine h1itr_prop_fcore(orb, torb)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : nbas, mval
  use mod_ormas, only : nfcore, nfun
  use mod_const, only : czero

  implicit none
  complex(c_double_complex), intent(in) :: orb(1:nbas, 1:nfun)
  complex(c_double_complex), intent(inout) :: torb(1:nbas, 1:nfun)

  complex(c_double_complex) :: tmp
  integer(c_long) :: ndfun, info, ifun, jfun
  integer(c_long), allocatable :: ipiv(:)
  complex(c_double_complex), allocatable :: amat(:,:)

 if (nfcore == 0) return

 allocate(ipiv(1:nfcore))
 allocate(amat(1:nfcore, 1:nfun))

 do ifun = 1, nfcore
    do jfun = 1, nfun
       if (mval(jfun) == mval(ifun)) then
          call zdotc_omp(nbas, orb(1, ifun), torb(1, jfun), tmp)
       else
          tmp = czero
       end if
       if (jfun <= nfcore) amat(ifun, jfun) = + tmp
       if (jfun >  nfcore) amat(ifun, jfun) = - tmp
    end do
 end do

 !debug write(6, "('amat-real:')")
 !debug do ifun = 1, nfcore
 !debug    do jfun = 1, nfun
 !debug       write(6, "(f12.8)", advance = 'no') dble(amat(ifun, jfun))
 !debug    end do
 !debug    write(6, *)
 !debug end do
 !debug write(6, "('amat-imag:')")
 !debug do ifun = 1, nfcore
 !debug    do jfun = 1, nfun
 !debug       write(6, "(f12.8)", advance = 'no') aimag(amat(ifun, jfun))
 !debug    end do
 !debug    write(6, *)
 !debug end do

 ndfun = nfun - nfcore
 call zgesv(nfcore, ndfun, amat, nfcore, ipiv, amat(1, nfcore + 1), nfcore, info)
 if (info /= 0) stop 'error in calling zgesv @ h1itr_prop.'

 do ifun = nfcore + 1, nfun
    do jfun = 1, nfcore
       if (mval(jfun) == mval(ifun)) then
!bug      call zaxpy_omp(nbas, amat(jfun, ifun), orb(1, jfun), torb(1, ifun))
          call zaxpy_omp(nbas, amat(jfun, ifun), torb(1, jfun), torb(1, ifun))
       end if
    end do
 end do

 deallocate(amat)
 deallocate(ipiv)

end subroutine h1itr_prop_fcore
!######################################################################
subroutine h1itr_prop_fcrot(lfield, dtime, wfn)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : nbas, mval
  use mod_ormas, only : nfcore, nfun
  use mod_const, only : half, iunit

  implicit none
  real(c_double), intent(in) :: lfield(9), dtime
  complex(c_double_complex), intent(inout) :: wfn(1:nbas, 1:nfun)
  complex(c_double_complex) :: zfac

 if (nfcore == 0) return

 zfac = exp(iunit * half * lfield(3) ** 2 * dtime)
 call zscal_omp(nbas * nfcore, zfac, wfn)

end subroutine h1itr_prop_fcrot
!!######################################################################
!subroutine h1itr_prop2(icomp, isplit, igauge, iprojfc, lfield, dtime, maxcyc, &
!                       thresh, cnpiv, cninv, orb, horb, torb, wfn)
!
!  use, intrinsic :: iso_c_binding
!  use mod_bas, only : nbas
!  use mod_sph, only : lmax1
!  use mod_rad, only : nrad, wrad
!  use mod_ormas, only : nfcore, nfun
!  use mod_const, only : zero, half, czero, runit, iunit
!  use mod_control, only : fedvr_normalized
!
!  implicit none
!  integer(c_long), intent(in) :: icomp, isplit, igauge, iprojfc, maxcyc, cnpiv(1:*)
!  real(c_double), intent(in) :: lfield(9), dtime, thresh
!  complex(c_double_complex), intent(in) :: cninv(1:*)
!  complex(c_double_complex), intent(out) :: orb(1:(nrad-1), 0:lmax1, 1:nfun)
!  complex(c_double_complex), intent(out) :: horb(1:(nrad-1), 0:lmax1, 1:nfun)
!  complex(c_double_complex), intent(out) :: torb(1:(nrad-1), 0:lmax1, 1:nfun)
!  complex(c_double_complex), intent(inout) :: wfn(1:(nrad-1), 0:lmax1, 1:nfun)
!
!  complex(c_double_complex) :: zfield, zfac1, zfac2, tmp
!  integer(c_long) :: ncyc, icyc, ifun, irad, l, size1, sizef, sized
!  logical, parameter :: debug = .true.
!  real(c_double) :: rmsres, maxres
!  logical projfc
!
!  projfc = iprojfc .ne. 0
!  size1 = nbas * nfun
!  sizef = nbas * nfcore
!  sized = size1 - sizef
!  zfield = lfield(3)
!  if (icomp == 1) then
!     zfac1 = - iunit * dtime * half
!     zfac2 = - iunit * dtime * half * zfield
!  else
!     zfac1 = - runit * dtime * half
!     zfac2 = - runit * dtime * half * zfield
!  end if
!
!  ! initialization
!  call zclear_omp(size1, horb)
!  call zcopy_omp(size1, wfn, orb)
!  call zcopy_omp(size1, wfn, torb)
!
!  !debug if (debug) then
!  !debug    write(6, "(i5, ':')", advance = 'no') -2
!  !debug    do ifun = 1, nfun
!  !debug       call zdotc_omp(nbas, orb(1, ifun), orb(1, ifun), tmp)
!  !debug       write(6, "(e20.5)", advance = 'no') sqrt(abs(tmp))
!  !debug    end do
!  !debug    write(6, *)
!  !debug end if
!
!  ! 1 - i * H * dt/2
!  call hprod_tprod_dyn(orb, horb)
!  if (isplit == 1 .and. igauge == 0) call hprod_zprod_dyn(zfield, orb, horb)
!  if (isplit == 1 .and. igauge == 1) call hprod_pzprod_dyn(zfield, orb, horb)
!
!  if (fedvr_normalized) then
!     call zaxpy_omp(sized, zfac1, horb(1,0,nfcore+1), orb(1,0,nfcore+1))
!  else
!     do ifun = nfcore + 1, nfun
!        do l = 0, lmax1
!           do irad = 1, nrad - 1
!              orb(irad, l, ifun) = orb(irad, l, ifun) + zfac1 * horb(irad, l, ifun) / wrad(irad)
!           end do
!        end do
!     end do
!  end if
!
!  !debug if (debug) then
!  !debug    write(6, "(i5, ':')", advance = 'no') -1
!  !debug    do ifun = 1, nfun
!  !debug       call zdotc_omp(nbas, orb(1, ifun), orb(1, ifun), tmp)
!  !debug       write(6, "(e20.5)", advance = 'no') sqrt(abs(tmp))
!  !debug    end do
!  !debug    write(6, *)
!  !debug end if
!
!  ! 1 / (1 + i * H0 dt /2)
!  call h1itr_cninv(cnpiv, cninv, orb)
!  call zcopy_omp(size1, orb, torb)
!
!  !debug if (debug) then
!  !debug    write(6, "(i5, ':')", advance = 'no') 0
!  !debug    do ifun = 1, nfun
!  !debug       call zdotc_omp(nbas, orb(1, ifun), orb(1, ifun), tmp)
!  !debug       write(6, "(e20.5)", advance = 'no') sqrt(abs(tmp))
!  !debug    end do
!  !debug    write(6, *)
!  !debug end if
!
!  if (isplit == 1) then
!     do icyc = 1, maxcyc
!        ncyc = icyc
!        maxres = zero
!        rmsres = zero
!        if (debug) write(6, "('h1itr_prop: ', i5, ':')", advance = 'no') icyc
!        do ifun = 1, nfun
!           if (fedvr_normalized) then
!              call zdotc_omp(nbas, orb(1,0,ifun), orb(1,0,ifun), tmp)
!           else
!              tmp = czero
!              do l = 0, lmax1
!                 do irad = 1, nrad - 1
!                    tmp = conjg(orb(irad, l, ifun)) * orb(irad, l, ifun) * wrad(irad)
!                 end do
!              end do
!           end if
!           if (debug) write(6, "(e20.5)", advance = 'no') sqrt(abs(tmp))
!
!           maxres = max(maxres, sqrt(abs(tmp)))
!           rmsres = rmsres + abs(tmp)
!        end do
!        if (debug) write(6, *)
!
!        rmsres = sqrt(rmsres / nfun)
!        if (rmsres < thresh) exit
!
!        ! i * V * dt/2
!        call zclear_omp(size1, horb)
!        if (projfc) then
!           if (igauge == 0) call hprod_zprod_all(zfac2, orb, horb);
!           if (igauge == 1) call hprod_pzprod_all(zfac2, orb, horb);
!        else
!           if (igauge == 0) call hprod_zprod_dyn(zfac2, orb, horb);
!           if (igauge == 1) call hprod_pzprod_dyn(zfac2, orb, horb);
!        end if
!
!        ! 1 / (1 + i * H0 dt /2)
!        if (fedvr_normalized) then
!           call zcopy_omp(size1, horb, orb);
!        else
!           do ifun = nfcore + 1, nfun
!              do l = 0, lmax1
!                 do irad = 1, nrad - 1
!                    orb(irad, l, ifun) = horb(irad, l, ifun) / wrad(irad)
!                 end do
!              end do
!           end do
!        end if
!
!        call h1itr_cninv(cnpiv, cninv, orb);
!        call zadd_omp(size1, orb, torb);
!
!     end do
!     if (debug) write(6, "('h1itr_prop: norm after ', i5, ' cycles: ' 2e20.5)") icyc, rmsres, maxres
!  end if
!
!  ! frozen-core correction
!  if (projfc) call h1itr_prop_fcore(wfn, torb)
!
!  call zcopy_omp(sized, torb(1,0,nfcore + 1), wfn(1,0,nfcore + 1))
!  !debug if (nfcore > 0 .and. igauge == 1) call h1itr_prop_fcrot(lfield, dtime, wfn)
!
!  !debug if (debug) then
!  !debug    write(6, "('orthonormalization check')")
!  !debug    do ifun = 1, nfun
!  !debug       do jfun = 1, nfun
!  !debug          if (mval(jfun) == mval(ifun)) then
!  !debug             call zdotc_omp(nbas, wfn(1, ifun), wfn(1, jfun), tmp)
!  !debug          else
!  !debug             tmp = czero
!  !debug          end if
!  !debug          write(6, "(2i5, 2e20.5)") ifun, jfun, tmp
!  !debug       end do
!  !debug    end do
!  !debug end if
!
!!  stop "@ h1itr_prop..."
!
!end subroutine h1itr_prop2
!######################################################################
