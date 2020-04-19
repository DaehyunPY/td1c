!######################################################################
subroutine dpade_prod(isplit, igauge, icomp, iprojfc, lfield, dtime, alpha, &
                      maxcyc, thresh, tpiv, tinv, wfn)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : nbas
  use mod_ormas, only : nfcore, nfun
  use mod_control, only : fedvr_normalized
  use mod_const, only : zero, half, czero, iunit, runit

  implicit none
  integer(c_long), intent(in) :: isplit, igauge, icomp, iprojfc, maxcyc, tpiv(1:*)
  real(c_double), intent(in) :: lfield(9), dtime, thresh
  complex(c_double_complex), intent(in) :: alpha
  complex(c_double_complex), intent(in) :: tinv(1:*)
  complex(c_double_complex), intent(inout) :: wfn(1:nbas, 1:nfun)

  complex(c_double_complex) :: zfac, zfac2, tmp
  integer(c_long) :: ncyc, icyc, ifun, size1, sized
  logical, parameter :: debug = .true.
  real(c_double) :: rmsres, maxres
  logical projfc
  complex(c_double_complex), allocatable :: orb(:,:)
  complex(c_double_complex), allocatable :: horb(:,:)
  complex(c_double_complex), allocatable :: torb(:,:)

  allocate(orb(1:nbas, 1:nfun))
  allocate(horb(1:nbas, 1:nfun))
  allocate(torb(1:nbas, 1:nfun))

  projfc = iprojfc == 1 .and. nfcore > 0
  size1 = nbas * nfun
  sized = size1 - nbas * nfcore
  if (icomp == 1) then
     zfac = - iunit * dtime / alpha * lfield(3)
     zfac2 = - iunit * dtime / alpha
  else
     zfac = - runit * dtime / alpha * lfield(3)
     zfac2 = - runit * dtime / alpha
  end if

  ! initialization
  call zclear_omp(size1, horb)
  call zcopy_omp(size1, wfn, orb)
  call zcopy_omp(size1, wfn, torb)

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

  if (isplit == 1) then
     do icyc = 1, maxcyc
        ncyc = icyc
        maxres = zero
        rmsres = zero
        if (debug) write(6, "('dpade_prod: cycle ', i5, ':')", advance = 'no') icyc
!       do ifun = 1, nfun
        do ifun = nfcore + 1, nfun
           call zdotc_omp(nbas, orb(1, ifun), orb(1, ifun), tmp)
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
     write(6, "('dpade_prod: norm after ', i5, ' cycles: ' 2e20.5)") icyc, rmsres, maxres
  end if

  ! frozen-core correction
  if (projfc) call dpade_prod_fcore(wfn, torb)

  call zcopy_omp(sized, torb(1, nfcore + 1), wfn(1, nfcore + 1))
  !debug if (nfcore > 0 .and. igauge == 1) call dpade_prod_fcrot(lfield, dtime, wfn)

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
!  stop "@ dpade_prod..."

end subroutine dpade_prod
!######################################################################
subroutine dpade_h0inv(tpiv, tinv, wfn)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad, ndvr
  use mod_sph, only : lmax1
  use mod_bas, only : mval
  use mod_ormas, only : nfun
  use mod_const, only : czero, runit

  implicit none
  integer(c_long), intent(in) :: tpiv(1:(nrad-1), 0:lmax1)
  complex(c_double_complex), intent(in) :: tinv(1:(3*ndvr+1), 1:(nrad-1), 0:lmax1)
  complex(c_double_complex), intent(inout) :: wfn(1:(nrad-1), 0:lmax1, 1:nfun)

  integer(c_long) :: ifun, l, dim, ld, info, lll, ull
  dim = nrad - 1
  ld = 3 * ndvr + 1

  !$omp parallel default(shared) private(lll, ull)
  !###########################
  call util_omp_disp(0, lmax1, lll, ull)
  do ifun = 1, nfun
     do l = max(lll, abs(mval(ifun))), ull
        call ZGBTRS('N', dim, ndvr, ndvr, 1, tinv(1,1,l), ld, tpiv(1,l), wfn(1,l,ifun), dim, info)
     end do
  end do
  !###########################
  !$omp end parallel

end subroutine dpade_h0inv
!######################################################################
subroutine dpade_prod_fcore(orb, torb)

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
 if (info /= 0) stop 'error in calling zgesv @ dpade_prod.'

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

end subroutine dpade_prod_fcore
!######################################################################
subroutine dpade_prod_fcrot(lfield, dtime, wfn)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : nbas
  use mod_ormas, only : nfcore, nfun
  use mod_const, only : half, iunit

  implicit none
  real(c_double), intent(in) :: lfield(3), dtime
  complex(c_double_complex), intent(inout) :: wfn(1:nbas, 1:nfun)
  complex(c_double_complex) :: zfac

 if (nfcore == 0) return

 zfac = exp(iunit * half * lfield(3) ** 2 * dtime)
 call zscal_omp(nbas * nfcore, zfac, wfn)

end subroutine dpade_prod_fcrot
!######################################################################
