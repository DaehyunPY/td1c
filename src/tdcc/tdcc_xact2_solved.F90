!################################################################################
subroutine tdcc_xact2_solved(dtime, cic, den1, bmat)
  !
  ! determine the time derivatives of ci coefficients and act-act rotations
  ! by direct inversion of coefficient matrix.
  ! bmat holds r.h.s on input, and solution on output.
  !
  use mod_control, only : icomp
  use mod_const, only : zero,czero,iunit
  use mod_ormas, only : iprint,nact,nrotaa,rotaa_mapb

  implicit none
  real(kind(0d0)), intent(in) :: dtime
  complex(kind(0d0)), intent(in) :: cic(1:*)
  complex(kind(0d0)), intent(in) :: den1(1:nact, 1:nact)
  complex(kind(0d0)), intent(inout) :: bmat(1:nact, 1:nact)
  !--------------------------------------------------------------------
  integer :: aact,iact,irot,numR,nrotR
  real(kind(0d0)), allocatable :: amat(:,:,:,:) ! coefficient matrix
  real(kind(0d0)), allocatable :: bvec(:,:)     ! r.h.s
  real(kind(0d0)), allocatable :: rvec(:,:)     ! solution

  if (icomp == 0) then
     numR = 1
  else
     numR = 2
  end if
  nrotR = numR*nrotaa

  allocate(bvec(1:nrotaa,1:nrotR))
  allocate(rvec(1:nrotaa,1:nrotR))
  allocate(amat(1:nrotaa,1:nrotR,1:nrotaa,1:nrotR))
  bvec(1:nrotaa,1:nrotR) = zero
  rvec(1:nrotaa,1:nrotR) = zero
  amat(1:nrotaa,1:nrotR,1:nrotaa,1:nrotR) = zero

  call tdcc_xact2_solved_bvec(numR, bmat, bvec)
  call tdcc_xact2_solved_amat(numR, cic, den1, amat)
  call tdcc_xact2_solved_solve(nrotR, amat, bvec, rvec)

  bmat(1:nact,1:nact) = czero
  do irot = 1, nrotaa
     aact = rotaa_mapb(irot,1)
     iact = rotaa_mapb(irot,2)

     if (iprint > 4) then
        if (icomp == 0) then
           write(6, "('bvec and rvec:',i5,2f20.10)") irot, bvec(irot,1), rvec(irot,1)
        else
           write(6, "('bvec and rvec:',i5,4f20.10)") irot, bvec(irot,1:2), rvec(irot,1:2)
        end if
     end if

     if (icomp == 0) then
        bmat(aact,iact) = rvec(irot,1)*dtime
        bmat(iact,aact) = -conjg(bmat(aact,iact))
     else
        bmat(aact,iact) = (rvec(irot,1)+iunit*rvec(irot,2))*dtime
        bmat(iact,aact) = -conjg(bmat(aact,iact))
     end if
  end do

  deallocate(rvec)
  deallocate(bvec)
  deallocate(amat)

end subroutine tdcc_xact2_solved
!################################################################################
!################################################################################
subroutine tdcc_xact2_solved_amat(numR, cic, den1, amat)

  use mod_control, only : icomp
  use mod_const, only : zero, czero, chalf
  use mod_cc, only : cc_rank
  use mod_ormas, only : iprint,nact,nrotaa,rotaa_mapb

  implicit none
  integer, intent(in) :: numR
  complex(kind(0d0)), intent(in) :: cic(1:*)
  complex(kind(0d0)), intent(in) :: den1(1:nact,1:nact)
  real(kind(0d0)), intent(out) :: amat(1:nrotaa,1:numR,1:nrotaa,1:numR)
  !--------------------------------------------------------------------
  integer :: aact,iact,bact,jact,irot,jrot
  complex(kind(0d0)) :: tmp
  complex(kind(0d0)), allocatable :: den2x(:,:,:,:)

  amat(1:nrotaa,1:numR,1:nrotaa,1:numR) = zero
  if (icomp == 0) then
     do irot = 1, nrotaa
        aact = rotaa_mapb(irot,1)
        iact = rotaa_mapb(irot,2)
        do jrot = 1, nrotaa
           bact = rotaa_mapb(jrot,1)
           jact = rotaa_mapb(jrot,2)
           if (aact == bact) amat(irot,1,jrot,1) = amat(irot,1,jrot,1) +  dble(den1(jact,iact))
           if (jact == iact) amat(irot,1,jrot,1) = amat(irot,1,jrot,1) -  dble(den1(aact,bact))
        end do
     end do
  else
     do irot = 1, nrotaa
        aact = rotaa_mapb(irot,1)
        iact = rotaa_mapb(irot,2)
        do jrot = 1, nrotaa
           bact = rotaa_mapb(jrot,1)
           jact = rotaa_mapb(jrot,2)
           if (aact == bact) then
              amat(irot,1,jrot,1) = amat(irot,1,jrot,1) +  dble(den1(jact,iact))
              amat(irot,1,jrot,2) = amat(irot,1,jrot,2) - aimag(den1(jact,iact))
              amat(irot,2,jrot,1) = amat(irot,2,jrot,1) + aimag(den1(jact,iact))
              amat(irot,2,jrot,2) = amat(irot,2,jrot,2) +  dble(den1(jact,iact))
           end if
           if (jact == iact) then
              amat(irot,1,jrot,1) = amat(irot,1,jrot,1) -  dble(den1(aact,bact))
              amat(irot,1,jrot,2) = amat(irot,1,jrot,2) + aimag(den1(aact,bact))
              amat(irot,2,jrot,1) = amat(irot,2,jrot,1) - aimag(den1(aact,bact))
              amat(irot,2,jrot,2) = amat(irot,2,jrot,2) -  dble(den1(aact,bact))
           end if
        end do
     end do
  end if

  if (cc_rank >= 3) then
     allocate(den2x(1:nact,1:nact,1:nact,1:nact))
     den2x = czero
     call tdcc_mkden2x(cic, den2x)
     if (icomp == 0) then
        do irot = 1, nrotaa
           aact = rotaa_mapb(irot,1)
           iact = rotaa_mapb(irot,2)
           do jrot = 1, nrotaa
              bact = rotaa_mapb(jrot,1)
              jact = rotaa_mapb(jrot,2)
              tmp = den2x(iact,aact,jact,bact) &
                  - den2x(jact,bact,iact,aact)
!##### I don't still know which is correct... #####
!             amat(irot,1,jrot,1) = amat(irot,1,jrot,1) - dble(tmp)
              amat(irot,1,jrot,1) = amat(irot,1,jrot,1) + dble(tmp)
!##### I don't still know which is correct... #####
           end do
        end do
     else
        do irot = 1, nrotaa
           aact = rotaa_mapb(irot,1)
           iact = rotaa_mapb(irot,2)
           do jrot = 1, nrotaa
              bact = rotaa_mapb(jrot,1)
              jact = rotaa_mapb(jrot,2)
              tmp = den2x(iact,aact,jact,bact) &
                  - den2x(jact,bact,iact,aact)
              amat(irot,1,jrot,1) = amat(irot,1,jrot,1) +  dble(tmp)
              amat(irot,1,jrot,2) = amat(irot,1,jrot,2) + aimag(tmp)
              amat(irot,2,jrot,1) = amat(irot,2,jrot,1) + aimag(tmp)
              amat(irot,2,jrot,2) = amat(irot,2,jrot,2) -  dble(tmp)
           end do
        end do
     end if
     deallocate(den2x)
  end if

  if (iprint > 4) then
     write(6, "('# ORMAS: amat11 ')")
     do irot = 1, nrotaa
        do jrot = 1, nrotaa
           write(6, "(e18.8)", advance='no') amat(irot,1,jrot,1)
        end do
        write(6, *)
     end do
     if (icomp == 1) then
        write(6, "('# ORMAS: amat12 ')")
        do irot = 1, nrotaa
           do jrot = 1, nrotaa
              write(6, "(e18.8)", advance='no') amat(irot,1,jrot,2)
           end do
           write(6, *)
        end do
        write(6, "('# ORMAS: amat21 ')")
        do irot = 1, nrotaa
           do jrot = 1, nrotaa
              write(6, "(e18.8)", advance='no') amat(irot,2,jrot,1)
           end do
           write(6, *)
        end do
        write(6, "('# ORMAS: amat22 ')")
        do irot = 1, nrotaa
           do jrot = 1, nrotaa
              write(6, "(e18.8)", advance='no') amat(irot,2,jrot,2)
           end do
           write(6, *)
        end do
     end if
  end if
  !stop 'for debug @ tdcc_xact2_solved_amat.'

end subroutine tdcc_xact2_solved_amat
!################################################################################
!################################################################################
subroutine tdcc_xact2_solved_bvec(numR, bmat, bvec)

  use mod_control, only : icomp
  use mod_const, only : zero, runit, iunit
  use mod_ormas, only : nact, nrotaa, rotaa_mapb

  implicit none
  integer, intent(in) :: numR
  complex(kind(0d0)), intent(in) :: bmat(1:nact, 1:nact)
  real(kind(0d0)), intent(out) :: bvec(1:nrotaa, 1:numR)
  !--------------------------------------------------------------------
  integer :: aact, iact, irot

  if (icomp == 0) then
     do irot = 1, nrotaa
        aact = rotaa_mapb(irot, 1)
        iact = rotaa_mapb(irot, 2)
        bvec(irot,1) = -dble(bmat(iact,aact))
     end do
  else
     do irot = 1, nrotaa
        aact = rotaa_mapb(irot, 1)
        iact = rotaa_mapb(irot, 2)
        bvec(irot,1) = aimag(bmat(iact,aact))
        bvec(irot,2) = -dble(bmat(iact,aact))
     end do
  end if

end subroutine tdcc_xact2_solved_bvec
!################################################################################
subroutine tdcc_xact2_solved_solve(nvar, amat, bvec, rvec)

  use mod_ormas, only : iprint
  use mod_control, only : reg_type, throcc3
  use mod_const, only : zero, one, two, three, half

  implicit none
  integer, intent(in) :: nvar
  real(kind(0d0)), intent(in) :: amat(1:nvar, 1:nvar)
  real(kind(0d0)), intent(in) :: bvec(1:nvar)
  real(kind(0d0)), intent(out) :: rvec(1:nvar)

  real(kind(0d0)) :: thresh
  real(kind(0d0)) :: inva, denom
  real(kind(0d0)), allocatable :: aeig(:)
  real(kind(0d0)), allocatable :: btmp(:)
  real(kind(0d0)), allocatable :: rtmp(:)
  real(kind(0d0)), allocatable :: uvec(:,:)
  real(kind(0d0)), allocatable :: vvec(:,:)
  integer :: ivar, jvar

  thresh = throcc3
  allocate(aeig(1:nvar))
  allocate(btmp(1:nvar))
  allocate(rtmp(1:nvar))
  allocate(uvec(1:nvar, 1:nvar))
  allocate(vvec(1:nvar, 1:nvar))

  call util_svd_real(nvar, nvar, amat, aeig, uvec, vvec)
  ! dpsi_smin = aeig(nvar)

  if (iprint > 4) then
     write(6, "('# tdcc_xact2_solved_solve: singular values of A')")
     do ivar = 1, min(8, nvar)
        write(6, "(i10, E15.5)") ivar, aeig(ivar)
     end do
  end if

  btmp(1:nvar) = zero
  do ivar = 1, nvar
     do jvar = 1, nvar
        btmp(ivar) = btmp(ivar) + uvec(jvar,ivar)*bvec(jvar)
     end do
  end do

  if (iprint > 4) then
     do ivar = 1, nvar
        write(6, "('b and U^t*b: ', i5, 4f20.10)") ivar, bvec(ivar), rvec(ivar)
     end do
  end if

  rtmp(1:nvar) = zero
  do ivar = 1, nvar
     if (reg_type == 0) then      ! 1/a <-- a/(a*a + d)
        denom = aeig(ivar) * aeig(ivar) + thresh
        inva = aeig(ivar) / denom
     else if (reg_type == 2) then ! 1/a <-- 1/(a + d * exp(-a / d))
        denom = aeig(ivar) + thresh * exp(-aeig(ivar) / thresh)
        inva = one / denom
!nyi     else if (reg_type == 3) then ! 1/a <-- 1/(a + d)
!nyi        denom = aeig(ivar) + thresh
!nyi        inva = one / denom
!nyi     else if (reg_type < 1) then  ! 1/a <-- 1/a
!nyi        denom = aeig(ivar)
!nyi        inva = one / denom
     else
        stop 'bad reg_type in tdcc_xact2_solved_solv.'
     end if

     aeig(ivar) = inva
     rtmp(ivar) = inva*btmp(ivar)
  end do

  rvec(1:nvar) = zero
  do ivar = 1, nvar
     do jvar = 1, nvar
        rvec(ivar) = rvec(ivar) + vvec(ivar, jvar) * rtmp(jvar)
     end do
  end do

!debug  do ivar = 1, nvar
!debug!     iact = rotaa_mapb(irot, 1)
!debug!     jact = rotaa_mapb(irot, 2)
!debug     write(6, "('bvec and rvec:: ', i5, 2f20.10)") ivar, bvec(ivar), rvec(ivar)
!debug  end do

  deallocate(vvec)
  deallocate(uvec)
  deallocate(rtmp)
  deallocate(btmp)
  deallocate(aeig)

end subroutine tdcc_xact2_solved_solve
!################################################################################
