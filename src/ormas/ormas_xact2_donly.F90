!################################################################################
subroutine ormas_xact2_donly(dtime, int1e, int2e, ene_act, cic, den1, den2, bmat, dcic)
  !
  ! determine the time derivatives of ci coefficients and act-act rotations
  ! by solving the coupled equations
  !
  use, intrinsic :: iso_c_binding
  use mod_control, only : icomp
  use mod_const, only : czero, runit, iunit
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
  integer(c_int) :: iact, jact
  complex(c_double_complex), allocatable :: bact(:,:) ! rhs: B = FD - DF*

  allocate(bact(1:nact, 1:nact))
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

     call ormas_xact2x_donly(dtime, cic, den1, den2, bact)

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

end subroutine ormas_xact2_donly
!################################################################################
subroutine ormas_xact2x_donly(dtime, cic, den1, den2, bmat)
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
  integer(c_int) :: iact, jact, irot
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

  call ormas_xact2x_donly_bvec(dtime, bmat, bvec, bnorm)
  call ormas_xact2x_donly_amat(cic, den1, den2, amat, anorm)
  if (icomp == 1) then
     call ormas_xact2x_solve(nrotaa, amat, bvec, rvec)
  else
     call ormas_xact2x_solve(nrotaa, amat, bvec, rvec)
     !call ormas_xact2x_donly_solver(nrotaa, amat, bvec, rvec)
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

end subroutine ormas_xact2x_donly
!################################################################################
subroutine ormas_xact2x_donly_amat(cic, den1, den2, amat, anorm)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : mval
  use mod_const, only : zero, one, czero
  use mod_ormas, only : ncore, nact, nrotaa, rotaa_mapb, ormas_nod2x, iprint

  implicit none
  complex(c_double_complex), intent(in) :: cic(1:*)
  complex(c_double_complex), intent(in) :: den1(1:nact, 1:nact)
  complex(c_double_complex), intent(in) :: den2(1:nact, 1:nact, 1:nact, 1:nact)
  complex(c_double), intent(out) :: amat(1:nrotaa, 1:nrotaa)
  real(c_double), intent(out) :: anorm
  !--------------------------------------------------------------------
  complex(c_double_complex), allocatable :: amat0(:,:,:,:)
  integer(c_int) :: iact, jact, kact, lact, irot, jrot
  complex(c_double_complex) :: tmp

  allocate(amat0(1:nact, 1:nact, 1:nact, 1:nact))
  amat0(1:nact, 1:nact, 1:nact, 1:nact) = czero

  do iact = 1, nact
     do jact = 1, nact
        if (mval(ncore+iact)/=mval(ncore+jact)) cycle
        do kact = 1, nact
           do lact = 1, nact
              if (mval(ncore+kact)/=mval(ncore+lact)) cycle
              if (iact == kact .and. mval(ncore+jact)==mval(ncore+lact)) &
                   amat0(iact, jact, kact, lact) &
               & = amat0(iact, jact, kact, lact) + den1(lact, jact)
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

  deallocate(amat0)

end subroutine ormas_xact2x_donly_amat
!################################################################################
subroutine ormas_xact2x_donly_bvec(dtime, bmat, bvec, bnorm)
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
  integer(c_int) :: iact, jact, irot
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
        write(6, "('ormas_xact2x_donly_bvec: ', 3i5, 2f20.10)") irot, iact, jact, bvec(irot)
     end if
     bnorm = bnorm + conjg(bvec(irot)) * bvec(irot)
  end do
  bnorm = sqrt(bnorm)

end subroutine ormas_xact2x_donly_bvec
!################################################################################
