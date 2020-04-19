!######################################################################
subroutine dpade_prod2(isplit, igauge, lfield, coeff1, maxcyc, thresh, tpiv, tinv, wfn)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : nbas
  use mod_ormas, only : nfcore, nfun
  use mod_control, only : fedvr_normalized
  use mod_hprod, only : orb, h1orb, torb
  use mod_const, only : zero, half, czero, iunit, runit

  implicit none
  real(c_double), intent(in) :: lfield(9), thresh
  complex(c_double_complex), intent(in) :: coeff1
  integer(c_int), intent(in) :: isplit, igauge, maxcyc, tpiv(1:*)
  complex(c_double_complex), intent(in) :: tinv(1:*)
  complex(c_double_complex), intent(inout) :: wfn(1:nbas, 1:nfun)

  complex(c_double_complex) :: zfac, zfac2, tmp
  integer(c_int) :: ncyc, icyc, ifun, size1, sized
  logical, parameter :: debug = .true.
  logical, parameter :: projfc = .false.
  real(c_double) :: rmsres, maxres

  size1 = nbas * nfun
  sized = size1 - nbas * nfcore
  zfac = -coeff1 * lfield(3)

  ! initialization
  call zclear_omp(size1, h1orb)
  call zcopy_omp(size1, wfn, orb)
  call zcopy_omp(size1, wfn, torb)

  ! 1 / (C0 + C1 * H0)
  call dpade_h0inv2(tpiv, tinv, orb)
  call zcopy_omp(size1, orb, torb)

  if (isplit == 1) then
     do icyc = 1, maxcyc
        ncyc = icyc
        maxres = zero
        rmsres = zero
        if (debug) write(6, "('dpade_prod: cycle ', i5, ':')", advance = 'no') icyc
        do ifun = nfcore + 1, nfun
           call zdotc_omp(nbas, orb(1, 0, ifun), orb(1, 0, ifun), tmp)
           if (debug) write(6, "(e20.5)", advance = 'no') sqrt(abs(tmp))

           maxres = max(maxres, sqrt(abs(tmp)))
           rmsres = rmsres + abs(tmp)
        end do
        if (debug) write(6, *)

        rmsres = sqrt(rmsres / nfun)
        if (rmsres < thresh) exit

        ! -C1 * V
        call zclear_omp(size1, h1orb)
        if (projfc) then
           if (igauge == 0) call hprod_zprod_all(zfac, orb, h1orb);
           if (igauge == 1) call hprod_pzprod_all(zfac, orb, h1orb);
           if (.not. fedvr_normalized) call hprod_invwrad_all(h1orb);
        else
           if (igauge == 0) call hprod_zprod_dyn(zfac, orb, h1orb);
           if (igauge == 1) call hprod_pzprod_dyn(zfac, orb, h1orb);
           if (.not. fedvr_normalized) call hprod_invwrad_dyn(h1orb);
        end if

        ! 1 / (C0 + C1 * H0)
        call zcopy_omp(size1, h1orb, orb);
        call dpade_h0inv2(tpiv, tinv, orb);
        call zadd_omp(size1, orb, torb);
     end do
     write(6, "('dpade_prod: norm after ', i5, ' cycles: ' 2e20.5)") icyc, rmsres, maxres
  end if

  ! frozen-core correction
  if (projfc) call dpade_prod_fcore2(wfn, torb)

  call zcopy_omp(sized, torb(1, 0, nfcore + 1), wfn(1, nfcore + 1))

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

!  stop "@ dpade_prod..."

end subroutine dpade_prod2
!######################################################################
subroutine dpade_h0inv2(tpiv, tinv, wfn)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad, ndvr
  use mod_sph, only : lmax1
  use mod_bas, only : mval
  use mod_ormas, only : nfun
  use mod_const, only : czero, runit

  implicit none
  integer(c_int), intent(in) :: tpiv(1:(nrad-1), 0:lmax1)
  complex(c_double_complex), intent(in) :: tinv(1:(3*ndvr+1), 1:(nrad-1), 0:lmax1)
  complex(c_double_complex), intent(inout) :: wfn(1:(nrad-1), 0:lmax1, 1:nfun)

  integer(c_int) :: ifun, l, dim, ld, info, lll, ull
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

end subroutine dpade_h0inv2
!######################################################################
subroutine dpade_prod_fcore2(orb, torb)

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
 if (info /= 0) stop 'error in calling zgesv @ dpade_prod.'

 do ifun = nfcore + 1, nfun
    do jfun = 1, nfcore
       if (mval(jfun) == mval(ifun)) then
          call zaxpy_omp(nbas, amat(jfun, ifun), torb(1, jfun), torb(1, ifun))
       end if
    end do
 end do

 deallocate(amat)
 deallocate(ipiv)

end subroutine dpade_prod_fcore2
!######################################################################
