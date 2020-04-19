!################################################################################
subroutine tdcc_xact2(dtime, int1e, int2e, den1, den2, cic, rmat, dcic)

  use, intrinsic :: iso_c_binding
  use mod_control, only : icomp, xact2_type, cionly
  use mod_const, only : one, czero, runit, iunit
  use mod_ormas, only : nfun,nact,ncore,ndetx
  use mod_cc, only : cc_rank,norb1,fock,int2x,optcc,bcc,optbcc

  implicit none
  real(kind(0d0)), intent(in) :: dtime
  complex(kind(0d0)), intent(in) :: int1e(1:nact, 1:nact)
  complex(kind(0d0)), intent(in) :: int2e(1:nact, 1:nact, 1:nact, 1:nact)
  complex(kind(0d0)), intent(in) :: den1(1:nact, 1:nact)
  complex(kind(0d0)), intent(in) :: den2(1:nact, 1:nact, 1:nact, 1:nact)
  complex(kind(0d0)), intent(in) :: cic(1:ndetx,1:2)
  complex(kind(0d0)), intent(inout) :: rmat(1:nfun, 1:nfun)
  complex(kind(0d0)), intent(inout) :: dcic(1:ndetx,1:2)
  !--------------------------------------------------------------------
  complex(kind(0d0)) :: fac1,fac2
  complex(kind(0d0)), allocatable :: xcic(:,:)
  complex(kind(0d0)), allocatable :: bmat(:,:)
  complex(kind(0d0)), allocatable :: bmat1(:,:)
  complex(kind(0d0)), allocatable :: bmat2(:,:)
  complex(kind(0d0)), allocatable :: bmat3(:,:)
  complex(kind(0d0)), allocatable :: soo(:,:), rsoo(:,:)
  complex(kind(0d0)), allocatable :: svv(:,:), rsvv(:,:)
  integer(c_int) :: h1,h2,h3,p1,p2,p3,iact,jact

  allocate(bmat(1:nact, 1:nact))
  bmat(1:nact, 1:nact) = czero

  call tdcc_mkint1x(int1e, int2e, fock)
  call tdcc_mkint2x(int2e, int2x)

  if (bcc .or. optbcc) then
     ! BCCD/BCCDT
     ! OBCCD/OBCCDT
     call tdcc_bcct1(cic, bmat)
     if (icomp == 1) then
        fock(:,:,1) = fock(:,:,1) - iunit*bmat(:,:)
        fock(:,:,2) = fock(:,:,2) - iunit*bmat(:,:)
     else
        fock(:,:,1) = fock(:,:,1) + bmat(:,:)
        fock(:,:,2) = fock(:,:,2) + bmat(:,:)
     end if
  end if

  call tdcc_hcc12(int1e,int2e,cic,dcic)

  !debug
  !write(6,"('tdcc_hcc12: hcc1')")
  !call tdcc_print(dcic)
  !call tdcc_hcc12(int1e,int2e,dcic,dcic)
  !write(6,"('tdcc_hcc12: hcc2')")
  !call tdcc_print(dcic)
  !stop
  !debug

  ! r.h.s of equations for a-a rotations
  if (optcc) then
     ! OCCD/OCCDT
     if (xact2_type == 0) then
        call tdcc_mkbmat(int1e,int2e,den1,den2,cic,dcic,bmat)
        call tdcc_xact2_solved(one, cic, den1, bmat)
     else 
        bmat = 0d0
        fac1 = +1d0
        fac2 = -1d0
        call tdcc_mkbmat1(fac1,int1e,int2e,den1,den2,bmat)
        call tdcc_mkbmat2(fac2,cic,dcic,bmat)
        call tdcc_xact2_solved_itr(one, cic, den1, bmat)
     end if
  else if (optbcc) then
     ! OBCCD/OBCCDT
     allocate(bmat1(1:nact,1:nact))
     allocate(bmat2(1:nact,1:nact))
     if (icomp == 1) then
        bmat1 = int1e -iunit*bmat
     else
        bmat1 = int1e + bmat
     end if
     bmat2 = 0d0
     call tdcc_mkbmat(bmat1,int2e,den1,den2,cic,dcic,bmat2)
     call tdcc_xact2_solveb(bmat2, cic, dcic)
     deallocate(bmat2)
     deallocate(bmat1)
  end if

  ! complete the orbital rotation matrix
  do iact = 1, nact
     do jact = 1, nact
        rmat(ncore+jact, ncore+iact) = bmat(jact, iact)*dtime
     end do
  end do

  ! real/imaginary time factor
  if (icomp == 1) then
     fac1 = -iunit * dtime
     fac2 =  iunit * dtime
  else
     fac1 = -runit * dtime
     fac2 = -runit * dtime
  end if
  dcic(1:ndetx,1) = dcic(1:ndetx,1) * fac1
  dcic(1:ndetx,2) = dcic(1:ndetx,2) * fac2

  ! a-a-rot contributions to ci derivatives, only for OCCDT
  if (cc_rank >= 3 .and. optcc .and. .not.cionly) then
     allocate(xcic(1:ndetx,1:2))
     call tdcc_hcc1(bmat, cic, xcic)
     if (icomp == 1) then
        dcic(1:ndetx,1) = dcic(1:ndetx,1) - xcic(1:ndetx,1)*dtime
        dcic(1:ndetx,2) = dcic(1:ndetx,2) + xcic(1:ndetx,2)*dtime
     else
        dcic(1:ndetx,1) = dcic(1:ndetx,1) - xcic(1:ndetx,1)*dtime
        dcic(1:ndetx,2) = dcic(1:ndetx,2) - xcic(1:ndetx,2)*dtime
     end if
     deallocate(xcic)
  end if

!  if (iprint > 4) then
!     write(6, "('# ORMAS: rmat-full')")
!     call util_matoutc(nfun, rmat)
!  end if

  deallocate(bmat)

end subroutine tdcc_xact2
!################################################################################
subroutine tdcc_xact2_test(dtime, int1e, int2e, den1, den2, cic, bmat, dcic)

  use, intrinsic :: iso_c_binding
  use mod_control, only : icomp
  use mod_const, only : czero, runit, iunit
  use mod_ormas, only : ormas_nocp, tdcc
  use mod_ormas, only : nfun, ncore, nact, nsub, ndetx, iprint
  use mod_cc, only : fock,int2x

  implicit none
  real(c_double), intent(in) :: dtime
  complex(c_double_complex), intent(in) :: int1e(1:*)
  complex(c_double_complex), intent(in) :: int2e(1:*)
  complex(c_double_complex), intent(in) :: den1(1:nact, 1:nact)
  complex(c_double_complex), intent(in) :: den2(1:nact, 1:nact, 1:nact, 1:nact)
  complex(c_double_complex), intent(in) :: cic(1:*)
  complex(c_double_complex), intent(inout) :: bmat(1:nfun,1:nfun)
  complex(c_double_complex), intent(inout) :: dcic(1:ndetx,1:2)
  !--------------------------------------------------------------------
  complex(c_double_complex) :: fac1,fac2
  complex(c_double_complex), allocatable :: dden(:,:) ! iunit times time-derivative of 1rdm
  complex(c_double_complex), allocatable :: bact(:,:) ! rhs: B = FD - DF* - dden(H)

  write(6,"('WARNING: tdcc_xact2ras: aa rotation neglected for tdcc.')")
  bmat = 0d0

  ! sigma vector
  call tdcc_hcc12(int1e, int2e, cic, dcic)

  if (icomp == 1) then
     fac1 = -iunit * dtime
     fac2 =  iunit * dtime
  else
     fac1 = -runit * dtime
     fac2 = -runit * dtime
  end if
  dcic(1:ndetx,1) = dcic(1:ndetx,1) * fac1
  dcic(1:ndetx,2) = dcic(1:ndetx,2) * fac2

end subroutine tdcc_xact2_test
!################################################################################
