!######################################################################
subroutine h1rat_prod(do_numer, do_denom, lfield, dtime, dimn, dimd, ncoeff, dcoeff1, tpiv, tinv, wfn)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : nbas
  use mod_sph, only : lmax1
  use mod_rad, only : nrad, ndvr
  use mod_ormas, only : nfcore, nfun
  use mod_hprod, only : orb, h1orb, torb
  use mod_const, only : zero, runit, cfalse
  use mod_control, only : icomp, isplit, iprojfc, h1rat_maxcyc, h1rat_thresh, PSP

  implicit none
  logical(c_bool), intent(in) :: do_numer, do_denom
  real(c_double), intent(in) :: lfield(9), dtime
  integer(c_int), intent(in) :: dimn, dimd
  complex(c_double_complex), intent(in) :: ncoeff(0:dimn)
  complex(c_double_complex), intent(in) :: dcoeff1(1:dimd)
  integer(c_int), intent(in) :: tpiv(1:(nrad-1), 0:lmax1, 1:dimd)
  complex(c_double_complex), intent(in) :: tinv(1:(3*ndvr+1), 1:(nrad-1), 0:lmax1, 1:dimd)
  complex(c_double_complex), intent(inout) :: wfn(1:(nrad-1), 0:lmax1, 1:*)

  complex(c_double_complex) :: zfac, tfac, tmp
  integer(c_int) :: idim, ncyc, icyc, ifun, idyn, ldyn
  logical, parameter :: debug = .false.
  real(c_double) :: rmsres, maxres
  logical projfc

  tfac = dtime
  idyn = nfcore + 1
  ldyn = nbas * (nfun - nfcore)
  projfc = iprojfc == 1 .and. nfcore > 0
  if (projfc) stop 'projfc not supported in h1rat_prod.'

  ! initialization
  call zcopy_omp(ldyn, wfn(1,0,idyn), orb(1,0,idyn))

  ! numerator
  if (do_numer) then
     call zcopy_omp(ldyn, wfn(1,0,idyn), torb(1,0,idyn))
     call zscal_omp(ldyn, ncoeff(0), orb(1,0,idyn))
     do idim = 1, dimn
        call zclear_omp(ldyn, h1orb(1,0,idyn))
        call hprod_h1tot(cfalse, lfield, torb, h1orb)
        call zscal2_omp(ldyn, tfac, h1orb(1,0,idyn), torb(1,0,idyn))
        call zaxpy_omp(ldyn, ncoeff(idim), torb(1,0,idyn), orb(1,0,idyn))
     end do
  end if

  ! denominator
  if (do_denom) then
     do idim = 1, dimd
        if (abs(dcoeff1(idim)) < h1rat_thresh) cycle
     
        ! 1 / (dcoeff0[i] + dcoeff1[i] * H0)
        !call h1rat_h0inv(tpiv(1,0,idim), tinv(1,1,0,idim), orb)
        call h1rat_h0inv_v2(tpiv(1,0,idim), tinv(1,1,0,idim), orb(1,0,1))
        ! matrix iteration method
        if (isplit == 1 .and. icomp == 1) then
           call zcopy_omp(ldyn, orb(1,0,idyn), torb(1,0,idyn))
           do icyc = 1, h1rat_maxcyc
              ! convergence check
              ncyc = icyc
              maxres = zero
              rmsres = zero
              if (debug) write(6, "('h1rat_prod: dim ', i5, ' cycle ', i5, ':')", advance = 'no') idim, icyc
              do ifun = nfcore + 1, nfun
                 call zdotc_omp(nbas, orb(1,0,ifun), orb(1,0,ifun), tmp)
                 if (debug) write(6, "(e20.5)", advance = 'no') sqrt(abs(tmp))
     
                 maxres = max(maxres, sqrt(abs(tmp)))
                 rmsres = rmsres + abs(tmp)
              end do
              if (debug) write(6, *)
     
              rmsres = sqrt(rmsres / nfun)
              if (rmsres < h1rat_thresh) exit
     
              ! ######################################################
              ! - dcoeff1[i] * V
              zfac = -dcoeff1(idim)
              call zclear_omp(ldyn, h1orb(1,0,idyn))
              call hprod_v1ext(cfalse, zfac, lfield, orb, h1orb);
              if (PSP) call hprod_projpp(zfac, lfield, orb, h1orb);
     
              ! 1 / (dcoeff0[i] + dcoeff1[i] * H0)
              call zcopy_omp(ldyn, h1orb(1,0,idyn), orb(1,0,idyn));
!              call h1rat_h0inv(tpiv(1,0,idim), tinv(1,1,0,idim), orb);
              call h1rat_h0inv_v2(tpiv(1,0,idim), tinv(1,1,0,idim), orb(1,0,1));
              call zadd_omp(ldyn, orb(1,0,idyn), torb(1,0,idyn));
              ! ######################################################
           end do
!           write(6, "('h1rat_prod: dim ', i5, ' convd ', i5, ':', 2e20.5)") idim, ncyc, rmsres, maxres
           call zcopy_omp(ldyn, torb(1,0,idyn), orb(1,0,idyn))
        end if
     end do
  end if

!  ! numerator
!  if (do_numer) then
!!     call zcopy_omp(ldyn, wfn(1,0,idyn), torb(1,0,idyn))
!     call zcopy_omp(ldyn, orb(1,0,idyn), torb(1,0,idyn))
!     call zscal_omp(ldyn, ncoeff(0), orb(1,0,idyn))
!     do idim = 1, dimn
!        call zclear_omp(ldyn, h1orb(1,0,idyn))
!        call hprod_h1tot(cfalse, lfield, torb, h1orb)
!        call zscal2_omp(ldyn, tfac, h1orb(1,0,idyn), torb(1,0,idyn))
!        call zaxpy_omp(ldyn, ncoeff(idim), torb(1,0,idyn), orb(1,0,idyn))
!     end do
!  end if

  call zcopy_omp(ldyn, orb(1,0,idyn), wfn(1,0,idyn))

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

!  stop "@ h1rat_prod..."

end subroutine h1rat_prod
!######################################################################
subroutine h1rat_prod_numer(lfield, dtime, dimn, zfactor, ncoeff, win, wout)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : nbas
  use mod_sph, only : lmax1
  use mod_rad, only : nrad, ndvr
  use mod_ormas, only : nfcore, nfun
  use mod_hprod, only : orb, h1orb, torb
  use mod_const, only : zero, runit, cfalse
  use mod_control, only : icomp, isplit, iprojfc, h1rat_maxcyc, h1rat_thresh

  implicit none
  real(c_double), intent(in) :: lfield(9), dtime
  integer(c_int), intent(in) :: dimn
  complex(c_double_complex), intent(in) :: zfactor
  complex(c_double_complex), intent(in) :: ncoeff(0:dimn)
  complex(c_double_complex), intent(in) :: win(1:(nrad-1), 0:lmax1, 1:*)
  complex(c_double_complex), intent(inout) :: wout(1:(nrad-1), 0:lmax1, 1:*)

  complex(c_double_complex) :: zfac, tfac, tmp
  integer(c_int) :: idim, ncyc, icyc, ifun, idyn, ldyn
  logical, parameter :: debug = .false.
  real(c_double) :: rmsres, maxres
  logical projfc

  tfac = dtime
  idyn = nfcore + 1
  ldyn = nbas * (nfun - nfcore)
  projfc = iprojfc == 1 .and. nfcore > 0
  if (projfc) stop 'projfc not supported in h1rat_prod.'

  ! initialization
  call zcopy_omp(ldyn, win(1,0,idyn), orb(1,0,idyn))

  ! numerator
  call zcopy_omp(ldyn, win(1,0,idyn), torb(1,0,idyn))
  call zscal_omp(ldyn, ncoeff(0), orb(1,0,idyn))
  do idim = 1, dimn
     call zclear_omp(ldyn, h1orb(1,0,idyn))
     call hprod_h1tot(cfalse, lfield, torb, h1orb)
     call zscal2_omp(ldyn, tfac, h1orb(1,0,idyn), torb(1,0,idyn))
     call zaxpy_omp(ldyn, ncoeff(idim), torb(1,0,idyn), orb(1,0,idyn))
  end do

  call zaxpy_omp(ldyn, zfactor, orb(1,0,idyn), wout(1,0,idyn))

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

!  stop "@ h1rat_prod..."

end subroutine h1rat_prod_numer
!######################################################################
subroutine h1rat_h0inv(tpiv, tinv, wfn)

  use, intrinsic :: iso_c_binding
  use mod_sph, only : lmax1
  use mod_bas, only : mval
  use mod_rad, only : nrad, ndvr
  use mod_ormas, only : nfun, nfcore, froz
  use mod_const, only : czero, runit

  implicit none
  integer(c_int), intent(in) :: tpiv(1:(nrad-1), 0:lmax1)
  complex(c_double_complex), intent(in) :: tinv(1:(3*ndvr+1), 1:(nrad-1), 0:lmax1)
  complex(c_double_complex), intent(inout) :: wfn(1:(nrad-1), 0:lmax1, 1:*)

  integer(c_int) :: ifun, l, dim, ld, info, lll, ull, irad
  dim = nrad - 1
  ld = 3 * ndvr + 1
  !$omp parallel default(shared) private(lll, ull)
  !###########################
  call util_omp_disp(0, lmax1, lll, ull)
  do ifun = nfcore + 1, nfun
     if (froz(ifun) < 0) cycle
     do l = max(lll, abs(mval(ifun))), ull
        call ZGBTRS('N', dim, ndvr, ndvr, 1, tinv(1,1,l), ld, tpiv(1,l), wfn(1,l,ifun), dim, info)
     end do
  end do
  !###########################
  !$omp end parallel

!  if (projhigh) then
!     call hprod_projhigh(wfn)
!  end if

end subroutine h1rat_h0inv
!######################################################################
subroutine h1rat_h0inv_v2(tpiv, tinv, wfn)

  use, intrinsic :: iso_c_binding
  use mod_sph, only : lmax1
  use mod_bas, only : mval
  use mod_rad, only : nrad, ndvr
  use mod_const, only : czero, runit
  use mod_ormas, only : nfcore, nfun, froz

  implicit none
  integer(c_int), intent(in) :: tpiv(1:(nrad-1), 0:lmax1)
  complex(c_double_complex), intent(in) :: tinv(1:(3*ndvr+1), 1:(nrad-1), 0:lmax1)
  complex(c_double_complex), intent(inout) :: wfn(1:(nrad-1), 0:lmax1, 1:*)

  integer(c_int) :: ifun, idyn, l, dim, ld, info, lll, ull, nproc, iproc, ndyn
  complex(c_double_complex), allocatable :: wtmp(:,:,:)
  integer(c_int), external :: util_omp_nproc
  integer(c_int), external :: util_omp_iproc
  dim = nrad - 1
  ld = 3 * ndvr + 1
  nproc = util_omp_nproc()
  ndyn = 0
  do ifun = nfcore + 1, nfun
     if (froz(ifun) >= 0) ndyn = ndyn + 1
  end do

  allocate(wtmp(1:(nrad-1), 1:ndyn, 0:(nproc-1)))
  
  !$omp parallel default(shared) private(iproc,idyn)
  iproc = util_omp_iproc()
  !$omp do
  !###########################
  do l = 0, lmax1
     idyn = 0
     do ifun = nfcore + 1, nfun
        if (froz(ifun) < 0) cycle
        idyn = idyn + 1
        wtmp(1:(nrad-1), idyn, iproc) = wfn(1:(nrad-1), l, ifun)
     end do

     call ZGBTRS('N', dim, ndvr, ndvr, ndyn, tinv(1,1,l), ld, tpiv(1,l), wtmp(1,1,iproc), dim, info)

     idyn = 0
     do ifun = nfcore + 1, nfun
        if (froz(ifun) < 0) cycle
        idyn = idyn + 1
        wfn(1:(nrad-1), l, ifun) = wtmp(1:(nrad-1), idyn, iproc)
     end do
  end do
  !###########################
  !$omp end do
  !$omp end parallel

  deallocate(wtmp)

!  if (projhigh) then
!     call hprod_projhigh(wfn)
!  end if

end subroutine h1rat_h0inv_v2
!######################################################################
