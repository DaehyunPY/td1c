!######################################################################
subroutine h1itr_prop2(icomp, isplit, igauge, lfield1, lfield2, dtime, maxcyc, &
     thresh, cnpiv, cninv, wfn)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : nbas
  use mod_sph, only : lmax1
  use mod_rad, only : nrad, wrad
  use mod_ormas, only : nfcore, nfun
  use mod_const, only : zero, half, czero, runit, iunit
  use mod_control, only : fedvr_normalized
  use mod_hprod, only : orb, orbg, torb, h1orb

  implicit none
  integer(c_int), intent(in) :: icomp, isplit, igauge, maxcyc, cnpiv(1:*)
  real(c_double), intent(in) :: lfield1(9), lfield2(9), dtime, thresh
  complex(c_double_complex), intent(in) :: cninv(1:*)
  complex(c_double_complex), intent(inout) :: wfn(1:(nrad-1), 0:lmax1, 1:nfun)

  complex(c_double_complex) :: zfield1, zfield2, tfac, zfac2, tmp
  integer(c_int) :: ncyc, icyc, ifun, size1, sizef, sized
  logical, parameter :: debug = .true.
  real(c_double) :: rmsres, maxres

  size1 = nbas * nfun
  sizef = nbas * nfcore
  sized = size1 - sizef
  zfield1 = lfield1(3)
  zfield2 = lfield2(3)
  if (icomp == 1) then
     tfac = - iunit * dtime * half
     zfac2 = - iunit * dtime * half * zfield2
  else
     tfac = - runit * dtime * half
     zfac2 = - runit * dtime * half * zfield2
  end if

  ! orbital at time t
  if (nfcore > 0) call hprod_fcorb(lfield1, orb, orbg)
  call zcopy_omp(sized, wfn(1,0,nfcore+1), orb(1,0,nfcore+1))

  ! 1 - i * Q_fc(t) * H * dt/2
  call zclear_omp(size1, h1orb)
  call hprod_tprod_dyn(orb, h1orb)
  if (isplit == 1 .and. igauge == 0) call hprod_zprod_dyn(zfield1, orb, h1orb)
  if (isplit == 1 .and. igauge == 1) call hprod_pzprod_dyn(zfield1, orb, h1orb)
  if (nfcore > 0) call hprod_projfc(.false., orb, h1orb)
  call zaxpy_omp(sized, tfac, h1orb(1,0,nfcore+1), orb(1,0,nfcore+1))

  ! 1 / (1 + i * H0 dt /2)
  if (nfcore > 0) call hprod_fcorb(lfield2, orb, orbg)
  call h1itr_cninv(cnpiv, cninv, orb)
  call zcopy_omp(size1, orb, torb)

  if (isplit == 1) then
     do icyc = 1, maxcyc
        ncyc = icyc
        maxres = zero
        rmsres = zero
        if (debug) write(6, "('h1itr_prop: ', i5, ':')", advance = 'no') icyc
        do ifun = 1, nfun
           call zdotc_omp(nbas, orb(1,0,ifun), orb(1,0,ifun), tmp)
           if (debug) write(6, "(e20.5)", advance = 'no') sqrt(abs(tmp))
           maxres = max(maxres, sqrt(abs(tmp)))
           rmsres = rmsres + abs(tmp)
        end do
        if (debug) write(6, *)

        rmsres = sqrt(rmsres / nfun)
        if (rmsres < thresh) exit

        ! i * V * dt/2
        call zclear_omp(size1, h1orb)
        if (igauge == 0) call hprod_zprod_all(zfac2, orb, h1orb);
        if (igauge == 1) call hprod_pzprod_all(zfac2, orb, h1orb);

        ! 1 / (1 + i * H0 dt /2)
        call zcopy_omp(size1, h1orb, orb);
        call h1itr_cninv(cnpiv, cninv, orb);
        call zadd_omp(size1, orb, torb);

     end do
     if (debug) write(6, "('h1itr_prop: norm after ', i5, ' cycles: ' 2e20.5)") icyc, rmsres, maxres
  end if

  ! frozen-core correction
  if (nfcore > 0) then
     call hprod_fcorb(lfield2, orb, orbg)
     call h1itr_prop2_fcore(orb, torb)
  end if

  call zcopy_omp(sized, torb(1,0,nfcore + 1), wfn(1,0,nfcore + 1))

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

!  stop "@ h1itr_prop..."

end subroutine h1itr_prop2
!######################################################################
subroutine h1itr_prop2_fcore(orb, torb)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : nbas, mval
  use mod_ormas, only : nfcore, nfun
  use mod_const, only : czero

  implicit none
  complex(c_double_complex), intent(in) :: orb(1:nbas, 1:nfun)
  complex(c_double_complex), intent(inout) :: torb(1:nbas, 1:nfun)

  complex(c_double_complex) :: tmp
  integer(c_int) :: ndfun, info, ifun, jfun
  integer(c_int), allocatable :: ipiv(:)
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

end subroutine h1itr_prop2_fcore
!######################################################################
