!######################################################################
subroutine dpade_cnic(isplit, igauge, icomp, iprojfc, lfield, dtime, alpha, &
                      maxcyc, thresh, tpiv, tinv, wfn)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : nbas
  use mod_sph, only : lmax1
  use mod_rad, only : nrad, wrad
  use mod_ormas, only : nfcore, nfun
  use mod_control, only : fedvr_normalized
  use mod_const, only : zero, half, czero, iunit, runit

  implicit none
  integer(c_long), intent(in) :: isplit, igauge, icomp, iprojfc, maxcyc, tpiv(1:*)
  real(c_double), intent(in) :: lfield(9), dtime, thresh
  complex(c_double_complex), intent(in) :: alpha
  complex(c_double_complex), intent(in) :: tinv(1:*)
  complex(c_double_complex), intent(inout) :: wfn(1:(nrad-1), 0:lmax1, 1:nfun)

  complex(c_double_complex) :: zfac, zfac2, zfield, tmp
  integer(c_long) :: ncyc, icyc, ifun, irad, l, size1, sized
  logical, parameter :: debug = .true.
  real(c_double) :: rmsres, maxres
  logical projfc
  complex(c_double_complex), allocatable :: orb(:,:,:)
  complex(c_double_complex), allocatable :: horb(:,:,:)
  complex(c_double_complex), allocatable :: torb(:,:,:)

  allocate(orb(1:(nrad-1), 0:lmax1, 1:nfun))
  allocate(horb(1:(nrad-1), 0:lmax1, 1:nfun))
  allocate(torb(1:(nrad-1), 0:lmax1, 1:nfun))

  projfc = iprojfc == 1 .and. nfcore > 0
  size1 = nbas * nfun
  sized = size1 - nbas * nfcore
  zfac = - iunit * dtime / alpha * lfield(3)
  zfac2 = - iunit * dtime / alpha
  zfield = lfield(3)

  ! initialization
  call zclear_omp(size1, horb)
  call zcopy_omp(size1, wfn, orb)
  call zcopy_omp(size1, wfn, torb)

  ! 1 - i * H * dt/2
  call hprod_tprod_dyn(orb, horb)
  if (isplit == 1 .and. igauge == 0) call hprod_zprod_dyn(zfield, orb, horb)
  if (isplit == 1 .and. igauge == 1) call hprod_pzprod_dyn(zfield, orb, horb)
  if (fedvr_normalized) then
     call zaxpy_omp(sized, zfac2, horb(1,0,nfcore+1), orb(1,0,nfcore+1))
  else
     do ifun = nfcore + 1, nfun
        do l = 0, lmax1
           do irad = 1, nrad - 1
              orb(irad, l, ifun) = orb(irad, l, ifun) + zfac2 * horb(irad, l, ifun) / wrad(irad)
           end do
        end do
     end do
  end if

  ! 1 / (1 + i * H0 dt /2)
  call dpade_h0inv(tpiv, tinv, orb)
  call zcopy_omp(size1, orb, torb)

  !debug if (debug) then
  !debug    write(6, "(i5, ':')", advance = 'no') 0
  !debug    do ifun = 1, nfun
  !debug       call zdotc_omp(nbas, orb(1, ifun), orb(1, ifun), tmp)
  !debug       write(6, "(e20.5)", advance = 'no') sqrt(abs(tmp))
  !debug    end do
  !debug    write(6, *)
  !debug end if

  if (icomp == 1 .and. isplit == 1) then
     do icyc = 1, maxcyc
        ncyc = icyc
        maxres = zero
        rmsres = zero
        if (debug) write(6, "('dpade_cnic: cycle ', i5, ':')", advance = 'no') icyc
!       do ifun = 1, nfun
        do ifun = nfcore + 1, nfun
           call zdotc_omp(nbas, orb(1,0,ifun), orb(1,0,ifun), tmp)
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
           if (igauge == 0) call hprod_zprod_all(zfac, orb, horb);
           if (igauge == 1) call hprod_pzprod_all(zfac, orb, horb);
           if (.not. fedvr_normalized) call hprod_invwrad_all(horb);
        else
           if (igauge == 0) call hprod_zprod_dyn(zfac, orb, horb);
           if (igauge == 1) call hprod_pzprod_dyn(zfac, orb, horb);
           if (.not. fedvr_normalized) call hprod_invwrad_dyn(horb);
        end if

        ! 1 / (1 + i * H0 dt /2)
        call zcopy_omp(size1, horb, orb);
        call dpade_h0inv(tpiv, tinv, orb);
        call zadd_omp(size1, orb, torb);
     end do
     write(6, "('dpade_cnic: norm after ', i5, ' cycles: ' 2e20.5)") icyc, rmsres, maxres
  end if

  ! frozen-core correction
  if (projfc) call dpade_prod_fcore(wfn, torb)

  call zcopy_omp(sized, torb(1,0,nfcore+1), wfn(1,0,nfcore+1))
  !debug if (nfcore > 0 .and. igauge == 1) call dpade_cnic_fcrot(lfield, dtime, wfn)

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
!  stop "@ dpade_cnic..."

end subroutine dpade_cnic
!######################################################################
