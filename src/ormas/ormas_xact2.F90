!################################################################################
subroutine ormas_xact2(dtime, int1e, int2e, ene_act, cic, den1, den2, bmat, dcic)
  !
  ! determine the time derivatives of ci coefficients and act-act rotations
  ! by solving the coupled equations
  !
  use, intrinsic :: iso_c_binding
  use mod_control, only : icomp
  use mod_const, only : czero, runit, iunit
  use mod_ormas, only : ormas_nocp, tdcc, ormas_donly
  use mod_ormas, only : nfun, ncore, nact, nsub, lcic, iprint

  implicit none
  real(c_double), intent(in) :: dtime
  complex(c_double_complex), intent(in) :: int1e(1:*)
  complex(c_double_complex), intent(in) :: int2e(1:*)
  real(c_double), intent(in) :: ene_act
  complex(c_double_complex), intent(in) :: cic(1:*)
  complex(c_double_complex), intent(in) :: den1(1:nact, 1:nact)
  complex(c_double_complex), intent(in) :: den2(1:nact, 1:nact, 1:nact, 1:nact)
  complex(c_double_complex), intent(inout) :: bmat(1:nfun, 1:nfun)
  complex(c_double_complex), intent(inout) :: dcic(1:*)
  !--------------------------------------------------------------------
  complex(c_double_complex) :: fac
  integer(c_long) :: iact, jact
  complex(c_double_complex), allocatable :: dden(:,:) ! iunit times time-derivative of 1rdm
  complex(c_double_complex), allocatable :: bact(:,:) ! rhs: B = FD - DF* - dden(H)

  if (tdcc) then
     call tdcc_xact2(dtime,int1e,int2e,den1,den2,cic,bmat,dcic)
     return
  else if (ormas_donly) then
     call ormas_xact2_donly(dtime, int1e, int2e, ene_act, cic, den1, den2, bmat, dcic)
     return
  end if

  allocate(dden(1:nact, 1:nact))
  allocate(bact(1:nact, 1:nact))

  dden(1:nact, 1:nact) = czero
  bact(1:nact, 1:nact) = czero

  ! sigma vector
  call ormas_hcic(int1e, int2e, cic, dcic, ene_act)

  if (nsub > 1) then
     ! r.h.s of equations for a-a rotations. I. FD - DF^{+}
     do iact = 1, nact
        do jact = 1, nact
           bact(jact, iact) = bmat(ncore+jact, ncore+iact)
        end do
     end do
     if (iprint > 4) then
        write(6, "('# ORMAS: bact-1')")
        call util_matoutc(nact, bact)
     end if

     ! r.h.s of equations for a-a rotations. II. -i dden(H)
     if (.not. ormas_nocp) call ormas_mkdden1(-runit, cic, dcic, bact)
     if (iprint > 4) then
        write(6, "('# ORMAS: bact-2')")
        call util_matoutc(nact, bact)
     end if

     call ormas_xact2x(dtime, cic, den1, den2, bact)

     if (iprint > 4) then
        write(6, "('# ORMAS: ract')")
        call util_matoutc(nact, bact)
     end if
  end if

  !########### Scale or clear <I|H|Psi> #############
  if (icomp == 1) then
     fac = -iunit * dtime
  else
     fac = -runit * dtime
  end if
  dcic(1:lcic) = dcic(1:lcic) * fac
  !##################################################

  if (nsub > 1) then
     ! a-a-rot contributions to ci derivatives
     fac = -runit
     call ormas_zhcic1(fac, bact, cic, dcic)

     ! complete the orbital rotation matrix
     do iact = 1, nact
        do jact = 1, nact
           bmat(ncore+jact, ncore+iact) = bact(jact, iact)
        end do
     end do
     if (iprint > 4) then
        write(6, "('# ORMAS: xmat-full')")
        call util_matoutc(nfun, bmat)
     end if
  end if

  deallocate(bact)
  deallocate(dden)

end subroutine ormas_xact2
!################################################################################
subroutine ormas_xact2x(dtime, cic, den1, den2, bmat)
  !
  ! determine the time derivatives of ci coefficients and act-act rotations
  ! by direct inversion of coefficient matrix.
  ! bmat holds r.h.s on input, and solution on output.
  !
  use, intrinsic :: iso_c_binding
  use mod_control, only : icomp
  use mod_const, only : zero, one, czero, runit, iunit
  use mod_ormas, only : ncore, nact, nrotaa, rotaa_mapb, iprint

  implicit none
  real(c_double), intent(in) :: dtime
  complex(c_double_complex), intent(in) :: cic(1:*)
  complex(c_double_complex), intent(in) :: den1(1:nact, 1:nact)
  complex(c_double_complex), intent(in) :: den2(1:nact, 1:nact, 1:nact, 1:nact)
  complex(c_double_complex), intent(inout) :: bmat(1:nact, 1:nact)
  !--------------------------------------------------------------------
  integer(c_long) :: iact, jact, irot
  real(c_double) :: bnorm, anorm, bmax
  complex(c_double_complex), allocatable :: amat(:,:) ! coefficient matrix
  complex(c_double_complex), allocatable :: bvec(:)   ! r.h.s
  complex(c_double_complex), allocatable :: rvec(:)   ! solution

  allocate(amat(1:nrotaa, 1:nrotaa))
  allocate(bvec(1:nrotaa))
  allocate(rvec(1:nrotaa))

  amat(1:nrotaa, 1:nrotaa) = czero
  bvec(1:nrotaa) = czero
  rvec(1:nrotaa) = czero

  call ormas_xact2x_bvec(dtime, bmat, bvec, bnorm)
  call ormas_xact2x_amat(cic, den1, den2, amat, anorm)
  if (icomp == 1) then
     call ormas_xact2x_solve(nrotaa, amat, bvec, rvec)
  else
     call ormas_xact2x_solve(nrotaa, amat, bvec, rvec)
     !call ormas_xact2x_solver(nrotaa, amat, bvec, rvec)
  end if

  bmax = zero
  bmat(1:nact, 1:nact) = czero
  do irot = 1, nrotaa
     iact = rotaa_mapb(irot, 1) ! greater
     jact = rotaa_mapb(irot, 2) ! lesser
     if (iprint > 4) then
        write(6, "('bvec and rvec:', i5, 4f20.10)") irot, bvec(irot), rvec(irot)
     end if
     bmat(iact, jact) =         rvec(irot)
     bmat(jact, iact) = - conjg(rvec(irot))
     bmax = max(bmax, abs(rvec(irot)))
  end do

  deallocate(rvec)
  deallocate(bvec)
  deallocate(amat)

end subroutine ormas_xact2x
!################################################################################
subroutine ormas_xact2x_amat(cic, den1, den2, amat, anorm)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : mval
  use mod_const, only : zero, one, czero
  use mod_ormas, only : ncore, nact, nrotaa, rotaa_mapb, ormas_nocp, ormas_nod2x, iprint

  implicit none
  complex(c_double_complex), intent(in) :: cic(1:*)
  complex(c_double_complex), intent(in) :: den1(1:nact, 1:nact)
  complex(c_double_complex), intent(in) :: den2(1:nact, 1:nact, 1:nact, 1:nact)
  complex(c_double), intent(out) :: amat(1:nrotaa, 1:nrotaa)
  real(c_double), intent(out) :: anorm
  !--------------------------------------------------------------------
  complex(c_double_complex), allocatable :: amat0(:,:,:,:)
  complex(c_double_complex), allocatable :: den2x(:,:,:,:)
  integer(c_long) :: iact, jact, kact, lact, irot, jrot
  complex(c_double_complex) :: tmp

  allocate(amat0(1:nact, 1:nact, 1:nact, 1:nact))
  allocate(den2x(1:nact, 1:nact, 1:nact, 1:nact))
  amat0(1:nact, 1:nact, 1:nact, 1:nact) = czero
  den2x(1:nact, 1:nact, 1:nact, 1:nact) = czero

  if (.not.(ormas_nocp .or. ormas_nod2x)) then
     call ormas_mkden2x(cic, den2x)
  end if

  do iact = 1, nact
     do jact = 1, nact
        if (mval(ncore+iact)/=mval(ncore+jact)) cycle
        do kact = 1, nact
           do lact = 1, nact
              if (mval(ncore+kact)/=mval(ncore+lact)) cycle
              if (iact == kact .and. mval(ncore+jact)==mval(ncore+lact)) &
                   amat0(iact, jact, kact, lact) &
               & = amat0(iact, jact, kact, lact) + den1(lact, jact)
              if (.not.(ormas_nocp .or. ormas_nod2x)) then
                 amat0(iact, jact, kact, lact) &
             & = amat0(iact, jact, kact, lact) &
             & + den2 (iact, jact, lact, kact) &
             & - den2x(iact, jact, kact, lact)
              end if
           end do
        end do
     end do
  end do

  do irot = 1, nrotaa
     iact = rotaa_mapb(irot, 1) ! greater
     jact = rotaa_mapb(irot, 2) ! lesser
     do jrot = 1, nrotaa
        kact = rotaa_mapb(jrot, 1) ! greater
        lact = rotaa_mapb(jrot, 2) ! lesser
        amat(irot, jrot) =  amat0(iact, jact, kact, lact)
     end do
  end do

  anorm = zero
  do irot = 1, nrotaa
     do jrot = 1, nrotaa
        tmp = amat(jrot, irot)
        anorm = anorm + dble(conjg(tmp) * tmp)
     end do
  end do
  anorm = sqrt(anorm)
        
  if (iprint > 4) then
     write(6, "('# ORMAS: amat-real')")
     do irot = 1, nrotaa
        do jrot = 1, nrotaa
           write(6, "(e12.4)", advance='no') dble(amat(jrot, irot))
        end do
        write(6, *)
     end do
     write(6, "('# ORMAS: amat-imag ')")
     do irot = 1, nrotaa
        do jrot = 1, nrotaa
           write(6, "(e12.4)", advance='no') aimag(amat(jrot, irot))
        end do
        write(6, *)
     end do
  end if

  deallocate(den2x)
  deallocate(amat0)

end subroutine ormas_xact2x_amat
!################################################################################
subroutine ormas_xact2x_bvec(dtime, bmat, bvec, bnorm)
  !
  use, intrinsic :: iso_c_binding
  use mod_control, only : icomp
  use mod_bas, only : mval
  use mod_const, only : zero, two, runit, iunit, ctwo, chalf
  use mod_ormas, only : ncore, nact, nrotaa, rotaa_mapb, iprint

  implicit none
  real(c_double), intent(in) :: dtime
  complex(c_double_complex), intent(in) :: bmat(1:nact, 1:nact)
  complex(c_double), intent(out) :: bvec(1:nrotaa)
  real(c_double), intent(out) :: bnorm
  !--------------------------------------------------------------------
  integer(c_long) :: iact, jact, irot
  complex(c_double_complex) :: tfac

  if (icomp == 1) then
     tfac = -iunit * dtime
  else
     tfac = -runit * dtime
  end if
!  tfac = runit

  bnorm = zero
  do irot = 1, nrotaa
     iact = rotaa_mapb(irot, 1) ! greater
     jact = rotaa_mapb(irot, 2) ! lesser
     if (mval(ncore+iact) == mval(ncore+jact)) then
        bvec(irot) = bmat(iact, jact) * tfac
     end if
     if (iprint > 4) then
        write(6, "('ormas_xact2x_bvec: ', 3i5, 2f20.10)") irot, iact, jact, bvec(irot)
     end if
     bnorm = bnorm + conjg(bvec(irot)) * bvec(irot)
  end do
  bnorm = sqrt(bnorm)

end subroutine ormas_xact2x_bvec
!################################################################################
subroutine ormas_xact2x_solve(nvar, amat, bvec, rvec)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : iprint
  use mod_const, only : czero, zero, one
  use mod_control, only : throcc3, reg_type, ncut_occ3

  implicit none
  integer(c_long), intent(in) :: nvar
  complex(c_double_complex), intent(inout) :: amat(1:nvar, 1:nvar)
  complex(c_double_complex), intent(in) :: bvec(1:nvar)
  complex(c_double_complex), intent(out) :: rvec(1:nvar)

  real(c_double) :: aeig, inva, denom, thresh
  complex(c_double_complex), allocatable :: btmp(:)
  complex(c_double_complex), allocatable :: rtmp(:)
  complex(c_double_complex), allocatable :: uvec(:,:)
  integer(c_long) :: ivar, jvar, ncutx

  thresh = throcc3
  allocate(btmp(1:nvar))
  allocate(rtmp(1:nvar))
  allocate(uvec(1:nvar, 1:nvar))

!debug  write(6, "('ormas_xact2x_solve: nvar = ', i10)") nvar
!debug  write(6, "('ormas_xact2x_solve: ncut = ', i10)") ormas_ncut

  call util_diag_comp(.true., nvar, amat, uvec)

  if (iprint > 1) then
     write(6, "('# Eigenvalues of A  ')", advance = 'no')
     do ivar = 1, nvar
        write(6, "(E15.5)", advance = 'no') dble(amat(ivar, ivar))
     end do
     write(6, *)
  end if

  btmp(1:nvar) = czero
  do ivar = 1, nvar
     do jvar = 1, nvar
        btmp(ivar) = btmp(ivar) + conjg(uvec(jvar, ivar)) * bvec(jvar)
     end do
  end do

  if (iprint > 4) then
     do ivar = 1, nvar
        write(6, "('b and U^+*b: ', i5, 4f20.10)") ivar, bvec(ivar), btmp(ivar)
     end do
  end if

  if (ncut_occ3 >= 0) then
     ncutx = ncut_occ3
  else
     ncutx = nvar
     write(6, "('# Eigenvalues of A  ')", advance = 'no')
     do ivar = 1, nvar
        write(6, "(E15.5)", advance = 'no') dble(amat(ivar, ivar))
     end do
     write(6, *)
  end if

  rtmp(1:nvar) = czero
  do ivar = 1, nvar - ncutx
!    if (abs(dble(amat(ivar, ivar))) < 1.E-15) then
     if (.false.) then
        amat(ivar, ivar) = zero
        rtmp(ivar) = zero
     else
        if (reg_type == 0) then ! 1/a <-- a/(a*a + d)
           aeig = dble(amat(ivar, ivar))
           denom = aeig * aeig + thresh
           inva = aeig / denom
!nyi     else if (reg_type == 1) then ! 1/a <-- a/(a*a + d)
!nyi        if (ivar <= nvar - ormas_ncut) then
!nyi           thresh = thresh1
!nyi        else
!nyi           thresh = thresh2
!nyi        end if
!nyi!        write(6, "('ivar = ', i5, ' thresh = ', E15.5)") ivar, thresh
!nyi        denom = aeig(ivar) * aeig(ivar) + thresh
!nyi        inva = aeig(ivar) / denom
        else if (reg_type == 2) then ! 1/a <-- 1/(a + d * exp(-a / d))
           aeig = dble(amat(ivar, ivar))
           denom = aeig + thresh * exp(-aeig / thresh)
           inva = one / denom
!nyi     else if (reg_type == 3) then ! 1/a <-- 1/(a + d)
!nyi        denom = aeig(ivar) + thresh2
!nyi        inva = one / denom
!nyi     else if (reg_type == 4) then ! 1/a <-- 1/max(d*max(a),a)
!nyi        if (ivar <= nvar - ormas_ncut) then
!nyi           denom = aeig(ivar)
!nyi           inva = one / denom
!nyi        else
!nyi           denom = max(thresh2 * aeig(1), aeig(ivar))
!nyi           inva = one / denom
!nyi        end if
!nyi     else if (reg_type == 5) then
!nyi        if (aeig(ivar) > min(thresh1, thresh2)) then
!nyi           inva = one / aeig(ivar)
!nyi!           write(6, "('ivar = ', i5, ' NO regularization: thresh = ', 2E15.5)") ivar, min(thresh1, thresh2), aeig(ivar)
!nyi        else
!nyi           inva = half * revthr * (three - (aeig(ivar) * revthr)**two)
!nyi!           write(6, "('ivar = ', i5, ' DO regularization: thresh = ', 2E15.5)") ivar, min(thresh1, thresh2), aeig(ivar)
!nyi        end if
!nyi     else if (reg_type < 0) then
!nyi        denom = aeig(ivar)
!nyi        inva = one / denom
        else
           stop 'bad reg_type in ormas_xact2x_solv.'
        end if
        amat(ivar, ivar) = one / inva
        rtmp(ivar) = inva * btmp(ivar)
     end if
  end do

!  if (iprint > 4) then
!     write(6, "('# Largest eigenvalues of mA ')", advance = 'no')
!     do ivar = 1, min(8, nvar)
!        write(6, "(E15.5)", advance = 'no') dble(amat(ivar, ivar))
!     end do
!     write(6, *)
!     write(6, "('# Smallest eigenvalues of mA')", advance = 'no')
!     do ivar = max(1, nvar - 7), nvar
!        write(6, "(E15.5)", advance = 'no') dble(amat(ivar, ivar))
!     end do
!     write(6, *)
!  end if

  rvec(1:nvar) = czero
  do ivar = 1, nvar
     do jvar = 1, nvar
        rvec(ivar) = rvec(ivar) + uvec(ivar, jvar) * rtmp(jvar)
     end do
  end do

  deallocate(uvec)
  deallocate(rtmp)
  deallocate(btmp)

!debug
!  stop 'for debug in ormas_xact2x_solve.'
!debug

end subroutine ormas_xact2x_solve
!################################################################################
subroutine ormas_xact2x_solver(nvar, amat, bvec, rvec)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : iprint
  use mod_const, only : zero, one
  use mod_control, only : throcc3, reg_type

  implicit none
  integer(c_long), intent(in) :: nvar
  complex(c_double_complex), intent(in) :: amat(1:nvar, 1:nvar)
  complex(c_double_complex), intent(in) :: bvec(1:nvar)
  complex(c_double_complex), intent(out) :: rvec(1:nvar)

  real(c_double) :: aeig, inva, denom, thresh
  real(c_double), allocatable :: atmp(:,:)
  real(c_double), allocatable :: btmp(:)
  real(c_double), allocatable :: rtmp(:)
  real(c_double), allocatable :: uvec(:,:)
  integer(c_long) :: ivar, jvar

  thresh = throcc3
  allocate(atmp(1:nvar, 1:nvar))
  allocate(btmp(1:nvar))
  allocate(rtmp(1:nvar))
  allocate(uvec(1:nvar, 1:nvar))

!debug  write(6, "('ormas_xact2x_solve: nvar = ', i10)") nvar
!debug  write(6, "('ormas_xact2x_solve: ncut = ', i10)") ormas_ncut

  atmp(1:nvar, 1:nvar) = dble(amat(1:nvar, 1:nvar))
  call util_diag_real(.true., nvar, atmp, uvec)

  if (iprint > 4) then
     write(6, "('# ormas_xact2x_solver: eigenvalues of A  ')", advance = 'no')
     do ivar = 1, nvar
        write(6, "(E15.5)", advance = 'no') atmp(ivar, ivar)
     end do
     write(6, *)
  end if

  btmp(1:nvar) = zero
  do ivar = 1, nvar
     do jvar = 1, nvar
        btmp(ivar) = btmp(ivar) + uvec(jvar, ivar) * dble(bvec(jvar))
     end do
  end do

  if (iprint > 4) then
     do ivar = 1, nvar
        write(6, "('b and U^+*b: ', i5, 4f20.10)") ivar, bvec(ivar), btmp(ivar)
     end do
  end if

  rtmp(1:nvar) = zero
  do ivar = 1, nvar
     if (reg_type == 0) then ! 1/a <-- a/(a*a + d)
        aeig = atmp(ivar, ivar)
        denom = aeig * aeig + thresh
        inva = aeig / denom
!nyi     else if (reg_type == 1) then ! 1/a <-- a/(a*a + d)
!nyi        if (ivar <= nvar - ormas_ncut) then
!nyi           thresh = thresh1
!nyi        else
!nyi           thresh = thresh2
!nyi        end if
!nyi!        write(6, "('ivar = ', i5, ' thresh = ', E15.5)") ivar, thresh
!nyi        denom = aeig(ivar) * aeig(ivar) + thresh
!nyi        inva = aeig(ivar) / denom
     else if (reg_type == 2) then ! 1/a <-- 1/(a + d * exp(-a / d))
        aeig = atmp(ivar, ivar)
        denom = aeig + thresh * exp(-aeig / thresh)
        inva = one / denom
!nyi     else if (reg_type == 3) then ! 1/a <-- 1/(a + d)
!nyi        denom = aeig(ivar) + thresh2
!nyi        inva = one / denom
!nyi     else if (reg_type == 4) then ! 1/a <-- 1/max(d*max(a),a)
!nyi        if (ivar <= nvar - ormas_ncut) then
!nyi           denom = aeig(ivar)
!nyi           inva = one / denom
!nyi        else
!nyi           denom = max(thresh2 * aeig(1), aeig(ivar))
!nyi           inva = one / denom
!nyi        end if
!nyi     else if (reg_type == 5) then
!nyi        if (aeig(ivar) > min(thresh1, thresh2)) then
!nyi           inva = one / aeig(ivar)
!nyi!           write(6, "('ivar = ', i5, ' NO regularization: thresh = ', 2E15.5)") ivar, min(thresh1, thresh2), aeig(ivar)
!nyi        else
!nyi           inva = half * revthr * (three - (aeig(ivar) * revthr)**two)
!nyi!           write(6, "('ivar = ', i5, ' DO regularization: thresh = ', 2E15.5)") ivar, min(thresh1, thresh2), aeig(ivar)
!nyi        end if
!nyi     else if (reg_type < 0) then
!nyi        denom = aeig(ivar)
!nyi        inva = one / denom
     else
        stop 'bad reg_type in ormas_xact2x_solv.'
     end if

     atmp(ivar, ivar) = one / inva
     rtmp(ivar) = inva * btmp(ivar)
  end do

  if (iprint > 4) then
     write(6, "('# Largest eigenvalues of mA ')", advance = 'no')
     do ivar = 1, min(8, nvar)
        write(6, "(E15.5)", advance = 'no') atmp(ivar, ivar)
     end do
     write(6, *)
     write(6, "('# Smallest eigenvalues of mA')", advance = 'no')
     do ivar = max(1, nvar - 7), nvar
        write(6, "(E15.5)", advance = 'no') atmp(ivar, ivar)
     end do
     write(6, *)
  end if

  rvec(1:nvar) = zero
  do ivar = 1, nvar
     do jvar = 1, nvar
        rvec(ivar) = rvec(ivar) + uvec(ivar, jvar) * rtmp(jvar)
     end do
  end do

  deallocate(uvec)
  deallocate(rtmp)
  deallocate(btmp)
  deallocate(atmp)

!debug
!  stop 'for debug in ormas_xact2x_solve.'
!debug

end subroutine ormas_xact2x_solver
!################################################################################
