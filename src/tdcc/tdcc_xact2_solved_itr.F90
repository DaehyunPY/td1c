!################################################################################
subroutine tdcc_xact2_solved_itr(dtime, cca, den1, bmat)
  !
  ! Determine the time derivatives of ci coefficients and act-act rotations
  ! by iterative inversion of coefficient matrix.
  ! bmat holds r.h.s without X contributions on input, and solution on output.
  !
  use, intrinsic :: iso_c_binding
  use mod_cc, only : cc_rank
  use mod_const, only : iunit,runit
  use mod_control, only : icomp,xact2_type,xact2_maxitr,xact2_thresh
  use mod_ormas, only : iprint,nact,nrotaa,rotaa_mapb,lcic
  use mod_cc, only : nXai,aiX,nXij,ijX,nXab,abX

  implicit none
  real(kind(0d0)), intent(in) :: dtime
  complex(kind(0d0)), intent(in) :: cca(1:*)
  complex(kind(0d0)), intent(in) :: den1(1:nact, 1:nact)
  complex(kind(0d0)), intent(inout) :: bmat(1:nact, 1:nact)
  !--------------------------------------------------------------------
  real(kind(0d0)) :: resd
  logical(c_bool) :: doitr
  integer(c_int) :: ai,p2,p1,irot,numR,nrotR,iitr
  complex(kind(0d0)),allocatable :: xcca(:)    ! -iX cont. i*dot(T)
  complex(kind(0d0)),allocatable :: bmath(:,:) ! H cont. to bmat
  complex(kind(0d0)),allocatable :: bmatx(:,:) ! -iX cont. to bmat
  real(kind(0d0)), allocatable :: bvec(:,:)     ! r.h.s
  real(kind(0d0)), allocatable :: rvec(:,:)     ! solution
  real(kind(0d0)),allocatable :: rvecP(:,:)    ! solution in the previous itr.
  real(kind(0d0)), allocatable :: ainv(:)
  real(kind(0d0)), allocatable :: uvec(:,:)
  real(kind(0d0)), allocatable :: vvec(:,:)
  real(kind(0d0)), allocatable :: amat(:,:,:,:)

  if (icomp == 0) then
     numR = 1
  else
     numR = 2
  end if
  nrotR = numR*nrotaa

  allocate(bmath(1:nact,1:nact))
  allocate(bmatx(1:nact,1:nact))
  allocate(xcca(1:lcic))
  allocate(bvec(1:nrotaa,1:numR))
  allocate(rvec(1:nrotaa,1:numR))
  allocate(rvecP(1:nrotaa,1:numR))
  allocate(ainv(1:nrotR))
  allocate(uvec(1:nrotR, 1:nrotR))
  allocate(vvec(1:nrotR, 1:nrotR))
  allocate(amat(1:nrotaa,1:numR,1:nrotaa,1:numR))

  amat = 0d0
  !call tdcc_xact2_solved_itr_amat(numR, cca, den1, amat)
  call tdcc_xact2_solved_itr_amat_general(numR, cca, den1, amat)
  call tdcc_xact2_solved_itr_ainv(nrotR, amat, ainv, uvec, vvec)

  bmath = bmat
  bmatx = 0d0
  rvecP = 0d0
  doitr = cc_rank>=3.and.xact2_type<=2

  do iitr = 0, xact2_maxitr
     bvec = 0d0
     rvec = 0d0
     call tdcc_xact2_solved_itr_bvec(numR, bmath, bmatx, bvec)
     call tdcc_xact2_solved_itr_solve(nrotR, ainv, uvec, vvec, bvec, rvec)

     resd = 0d0
     do irot = 1,nrotaa
        resd = resd + dble((rvec(irot,1)-rvecP(irot,1))**2)
        if (icomp == 1) resd = resd + dble((rvec(irot,2)-rvecP(irot,2))**2)
     end do
     resd = sqrt(abs(resd))
     if (iprint > 1) then
        write(6,"('tdcc_xact2_solved_itr:',i5,f20.10)") iitr,resd
     end if

     if (abs(resd) < xact2_thresh) exit
     rvecP = rvec

     bmat = 0d0
     do irot = 1, nrotaa
        p2 = rotaa_mapb(irot,1)
        p1 = rotaa_mapb(irot,2)
        if (iprint > 4) then
           if (icomp == 0) then
              write(6, "('bvec and rvec:',3i5,2f20.10)") &
                   irot,p1,p2,bvec(irot,1),rvec(irot,1)
           else
              write(6, "('bvec and rvec:',3i5,4f20.10)") &
                   irot,p1,p2,bvec(irot,1:2),rvec(irot,1:2)
           end if
        end if
        if (icomp == 0) then
           bmat(p2,p1) = rvec(irot,1)
           bmat(p1,p2) = -conjg(bmat(p2,p1))
        else
           bmat(p2,p1) = (rvec(irot,1)+iunit*rvec(irot,2))
           bmat(p1,p2) = -conjg(bmat(p2,p1))
        end if
     end do

     if (doitr) then
        xcca = 0d0;
        bmatx = bmat;
        if (xact2_type == 2) then
           do ai = 1, nXai
              irot = aiX(ai)
              p2 = rotaa_mapb(irot,1)
              p1 = rotaa_mapb(irot,2)
              bmatx(p2,p1) = 0d0
              bmatx(p1,p2) = 0d0
           end do
        end if
        call tdcc_hcc1(bmatx,cca,xcca)

        bmatx = 0d0;
        call tdcc_mkbmat2(runit,cca,xcca,bmatx)
     else
        exit
     end if
  end do

  bmat = bmat*dtime

  deallocate(amat)
  deallocate(ainv)
  deallocate(uvec)
  deallocate(vvec)
  deallocate(rvecP)
  deallocate(rvec)
  deallocate(bvec)
  deallocate(xcca)
  deallocate(bmatx)
  deallocate(bmath)

end subroutine tdcc_xact2_solved_itr
!################################################################################
!################################################################################
subroutine tdcc_xact2_solved_itr_amat(numR, cca, den1, amat)

  use, intrinsic :: iso_c_binding
  use mod_control, only : xact2_type
  use mod_cc, only : cc_rank,nXai,aiX,nXij,ijX,nXab,abX
  use mod_ormas, only : iprint,nact,nrotaa,rotaa_mapb

  implicit none
  integer(c_int), intent(in) :: numR
  complex(kind(0d0)), intent(in) :: cca(1:*)
  complex(kind(0d0)), intent(in) :: den1(1:nact,1:nact)
  real(kind(0d0)), intent(inout) :: amat(1:nrotaa,1:numR,1:nrotaa,1:numR)
  !--------------------------------------------------------------------
  integer(c_int) :: i,j,a,b,ai,bj,i1,i2,j1,j2,a1,a2,b1,b2,i1i2,j1j2,a1a2,b1b2,irot,jrot
  complex(kind(0d0)) :: Aval,Bval
  complex(kind(0d0)), allocatable :: AB(:,:,:,:)

  if (cc_rank >= 3 .and. xact2_type >= 2) then
     allocate(AB(1:nact,1:nact,1:nact,1:nact))
     AB = 0d0
     call tdcc_mkamat(cca, AB)
     if (xact2_type == 2) then
        do ai = 1, nXai
           irot = aiX(ai)
           a = rotaa_mapb(irot,1)
           i = rotaa_mapb(irot,2)
           do bj = 1, nXai
              jrot = aiX(bj)
              b = rotaa_mapb(jrot,1)
              j = rotaa_mapb(jrot,2)
              Bval = AB(a,i,b,j) - AB(b,j,a,i)
              amat(irot,1,jrot,1) = amat(irot,1,jrot,1) +  dble(Bval)
              if (numR == 2) amat(irot,1,jrot,2) = amat(irot,1,jrot,2) + aimag(Bval)
              if (numR == 2) amat(irot,2,jrot,1) = amat(irot,2,jrot,1) + aimag(Bval)
              if (numR == 2) amat(irot,2,jrot,2) = amat(irot,2,jrot,2) -  dble(Bval)
           end do
        end do
     end if
     deallocate(AB)
  end if

  do ai = 1, nXai
     irot = aiX(ai)
     a = rotaa_mapb(irot,1)
     i = rotaa_mapb(irot,2)
     do bj = 1, nXai
        jrot = aiX(bj)
        b = rotaa_mapb(jrot,1)
        j = rotaa_mapb(jrot,2)
        if (a == b) then
           Aval = +den1(j,i)
           amat(irot,1,jrot,1) = amat(irot,1,jrot,1) +  dble(Aval)
           if (numR == 2) amat(irot,1,jrot,2) = amat(irot,1,jrot,2) - aimag(Aval)
           if (numR == 2) amat(irot,2,jrot,1) = amat(irot,2,jrot,1) + aimag(Aval)
           if (numR == 2) amat(irot,2,jrot,2) = amat(irot,2,jrot,2) +  dble(Aval)
        end if
        if (i == j) then
           Aval = -den1(a,b)
           amat(irot,1,jrot,1) = amat(irot,1,jrot,1) +  dble(Aval)
           if (numR == 2) amat(irot,1,jrot,2) = amat(irot,1,jrot,2) - aimag(Aval)
           if (numR == 2) amat(irot,2,jrot,1) = amat(irot,2,jrot,1) + aimag(Aval)
           if (numR == 2) amat(irot,2,jrot,2) = amat(irot,2,jrot,2) +  dble(Aval)
        end if
     end do
     do j1j2 = 1, nXij
        jrot = ijX(j1j2)
        j1 = rotaa_mapb(jrot,1)
        j2 = rotaa_mapb(jrot,2)
        if (i == j2) then
           Aval = -den1(a,j1)
           amat(irot,1,jrot,1) = amat(irot,1,jrot,1) +  dble(Aval)
           if (numR == 2) amat(irot,1,jrot,2) = amat(irot,1,jrot,2) - aimag(Aval)
           if (numR == 2) amat(irot,2,jrot,1) = amat(irot,2,jrot,1) + aimag(Aval)
           if (numR == 2) amat(irot,2,jrot,2) = amat(irot,2,jrot,2) +  dble(Aval)
        end if
        if (i == j1) then
           Bval = +den1(a,j2)
           amat(irot,1,jrot,1) = amat(irot,1,jrot,1) +  dble(Bval)
           if (numR == 2) amat(irot,1,jrot,2) = amat(irot,1,jrot,2) + aimag(Bval)
           if (numR == 2) amat(irot,2,jrot,1) = amat(irot,2,jrot,1) + aimag(Bval)
           if (numR == 2) amat(irot,2,jrot,2) = amat(irot,2,jrot,2) -  dble(Bval)
        end if
     end do
     do b1b2 = 1, nXab
        jrot = abX(b1b2)
        b2 = rotaa_mapb(jrot,1)
        b1 = rotaa_mapb(jrot,2)
        if (a == b1) then
           Bval = -den1(b2,i)
           amat(irot,1,jrot,1) = amat(irot,1,jrot,1) +  dble(Bval)
           if (numR == 2) amat(irot,1,jrot,2) = amat(irot,1,jrot,2) + aimag(Bval)
           if (numR == 2) amat(irot,2,jrot,1) = amat(irot,2,jrot,1) + aimag(Bval)
           if (numR == 2) amat(irot,2,jrot,2) = amat(irot,2,jrot,2) -  dble(Bval)
        end if
        if (a == b2) then
           Aval = +den1(b1,i)
           amat(irot,1,jrot,1) = amat(irot,1,jrot,1) +  dble(Aval)
           if (numR == 2) amat(irot,1,jrot,2) = amat(irot,1,jrot,2) - aimag(Aval)
           if (numR == 2) amat(irot,2,jrot,1) = amat(irot,2,jrot,1) + aimag(Aval)
           if (numR == 2) amat(irot,2,jrot,2) = amat(irot,2,jrot,2) +  dble(Aval)
        end if
     end do
  end do
  do i1i2 = 1, nXij
     irot = ijX(i1i2)
     i1 = rotaa_mapb(irot,1)
     i2 = rotaa_mapb(irot,2)
     do bj = 1, nXai
        jrot = aiX(bj)
        b = rotaa_mapb(jrot,1)
        j = rotaa_mapb(jrot,2)
        if (i1 == j) then
           Bval = -den1(b,i2)
           amat(irot,1,jrot,1) = amat(irot,1,jrot,1) +  dble(Bval)
           if (numR == 2) amat(irot,1,jrot,2) = amat(irot,1,jrot,2) + aimag(Bval)
           if (numR == 2) amat(irot,2,jrot,1) = amat(irot,2,jrot,1) + aimag(Bval)
           if (numR == 2) amat(irot,2,jrot,2) = amat(irot,2,jrot,2) -  dble(Bval)
        end if
        if (i2 == j) then
           Aval = -den1(i1,b)
           amat(irot,1,jrot,1) = amat(irot,1,jrot,1) +  dble(Aval)
           if (numR == 2) amat(irot,1,jrot,2) = amat(irot,1,jrot,2) - aimag(Aval)
           if (numR == 2) amat(irot,2,jrot,1) = amat(irot,2,jrot,1) + aimag(Aval)
           if (numR == 2) amat(irot,2,jrot,2) = amat(irot,2,jrot,2) +  dble(Aval)
        end if
     end do
     do j1j2 = 1, nXij
        jrot = ijX(j1j2)
        j1 = rotaa_mapb(jrot,1)
        j2 = rotaa_mapb(jrot,2)
        if (i1 == j1) then
           Aval = +den1(j2,i2)
           amat(irot,1,jrot,1) = amat(irot,1,jrot,1) +  dble(Aval)
           if (numR == 2) amat(irot,1,jrot,2) = amat(irot,1,jrot,2) - aimag(Aval)
           if (numR == 2) amat(irot,2,jrot,1) = amat(irot,2,jrot,1) + aimag(Aval)
           if (numR == 2) amat(irot,2,jrot,2) = amat(irot,2,jrot,2) +  dble(Aval)
        end if
        if (i2 == j2) then
           Aval = -den1(i1,j1)
           amat(irot,1,jrot,1) = amat(irot,1,jrot,1) +  dble(Aval)
           if (numR == 2) amat(irot,1,jrot,2) = amat(irot,1,jrot,2) - aimag(Aval)
           if (numR == 2) amat(irot,2,jrot,1) = amat(irot,2,jrot,1) + aimag(Aval)
           if (numR == 2) amat(irot,2,jrot,2) = amat(irot,2,jrot,2) +  dble(Aval)
        end if
     end do
  end do
  do a1a2 = 1, nXab
     irot = abX(a1a2)
     a2 = rotaa_mapb(irot,1)
     a1 = rotaa_mapb(irot,2)
     do bj = 1, nXai
        jrot = aiX(bj)
        b = rotaa_mapb(jrot,1)
        j = rotaa_mapb(jrot,2)
        if (a2 == b) then
           Aval = +den1(j,a1)
           amat(irot,1,jrot,1) = amat(irot,1,jrot,1) +  dble(Aval)
           if (numR == 2) amat(irot,1,jrot,2) = amat(irot,1,jrot,2) - aimag(Aval)
           if (numR == 2) amat(irot,2,jrot,1) = amat(irot,2,jrot,1) + aimag(Aval)
           if (numR == 2) amat(irot,2,jrot,2) = amat(irot,2,jrot,2) +  dble(Aval)
        end if
        if (a1 == b) then
           Bval = +den1(a2,j)
           amat(irot,1,jrot,1) = amat(irot,1,jrot,1) +  dble(Bval)
           if (numR == 2) amat(irot,1,jrot,2) = amat(irot,1,jrot,2) + aimag(Bval)
           if (numR == 2) amat(irot,2,jrot,1) = amat(irot,2,jrot,1) + aimag(Bval)
           if (numR == 2) amat(irot,2,jrot,2) = amat(irot,2,jrot,2) -  dble(Bval)
        end if
     end do
     do b1b2 = 1, nXab
        jrot = abX(b1b2)
        b2 = rotaa_mapb(jrot,1)
        b1 = rotaa_mapb(jrot,2)
        if (a2 == b2) then
           Aval = +den1(b1,a1)
           amat(irot,1,jrot,1) = amat(irot,1,jrot,1) +  dble(Aval)
           if (numR == 2) amat(irot,1,jrot,2) = amat(irot,1,jrot,2) - aimag(Aval)
           if (numR == 2) amat(irot,2,jrot,1) = amat(irot,2,jrot,1) + aimag(Aval)
           if (numR == 2) amat(irot,2,jrot,2) = amat(irot,2,jrot,2) +  dble(Aval)
        end if
        if (a1 == b1) then
           Aval = -den1(a2,b2)
           amat(irot,1,jrot,1) = amat(irot,1,jrot,1) +  dble(Aval)
           if (numR == 2) amat(irot,1,jrot,2) = amat(irot,1,jrot,2) - aimag(Aval)
           if (numR == 2) amat(irot,2,jrot,1) = amat(irot,2,jrot,1) + aimag(Aval)
           if (numR == 2) amat(irot,2,jrot,2) = amat(irot,2,jrot,2) +  dble(Aval)
        end if
     end do
  end do

end subroutine tdcc_xact2_solved_itr_amat
!################################################################################
subroutine tdcc_xact2_solved_itr_amat_general(numR, cca, den1, amat)

  use, intrinsic :: iso_c_binding
  use mod_control, only : xact2_type
  use mod_cc, only : cc_rank,nXai,aiX,nXij,ijX,nXab,abX
  use mod_ormas, only : iprint,nact,nrotaa,rotaa_mapb

  implicit none
  integer(c_int), intent(in) :: numR
  complex(kind(0d0)), intent(in) :: cca(1:*)
  complex(kind(0d0)), intent(in) :: den1(1:nact,1:nact)
  real(kind(0d0)), intent(inout) :: amat(1:nrotaa,1:numR,1:nrotaa,1:numR)
  !--------------------------------------------------------------------
  integer(c_int) :: irot,jrot,p1,p2,q1,q2
  real(kind(0d0)) :: AR,AI,BR,BI
  complex(kind(0d0)), allocatable :: alpha(:,:,:,:)
  complex(kind(0d0)), allocatable :: cmat(:,:,:,:)

  allocate(cmat(1:nact,1:nact,1:nact,1:nact))
  cmat = 0d0

  do irot = 1,nrotaa
     p2 = rotaa_mapb(irot,1)
     p1 = rotaa_mapb(irot,2)
     do jrot = 1,nrotaa
        q2 = rotaa_mapb(jrot,1)
        q1 = rotaa_mapb(jrot,2)
        if (p2 == q2) cmat(p2,p1,q1,q2) = cmat(p2,p1,q1,q2) + den1(q1,p1)
        if (p1 == q1) cmat(p2,p1,q1,q2) = cmat(p2,p1,q1,q2) - den1(p2,q2)
        if (p2 == q1) cmat(p2,p1,q2,q1) = cmat(p2,p1,q2,q1) + den1(q2,p1)
        if (p1 == q2) cmat(p2,p1,q2,q1) = cmat(p2,p1,q2,q1) - den1(p2,q1)
     end do
  end do

  if (cc_rank >= 3 .and. xact2_type >= 2) then
     allocate(alpha(1:nact,1:nact,1:nact,1:nact))
     alpha = 0d0
     call tdcc_mkamat(cca, alpha)
     do irot = 1,nrotaa
        p2 = rotaa_mapb(irot,1)
        p1 = rotaa_mapb(irot,2)
        do jrot = 1,nrotaa
           q2 = rotaa_mapb(jrot,1)
           q1 = rotaa_mapb(jrot,2)
           cmat(p2,p1,q1,q2) = cmat(p2,p1,q1,q2) &
                -(alpha(p2,p1,q1,q2) - alpha(q1,q2,p2,p1) &
          - conjg(alpha(p1,p2,q2,q1) - alpha(q2,q1,p1,p2)))
           cmat(p2,p1,q2,q1) = cmat(p2,p1,q2,q1) &
                -(alpha(p2,p1,q2,q1) - alpha(q2,q1,p2,p1) &
          - conjg(alpha(p1,p2,q1,q2) - alpha(q1,q2,p1,p2)))
        end do
     end do
     deallocate(alpha)
  end if

  do irot = 1,nrotaa
     p2 = rotaa_mapb(irot,1)
     p1 = rotaa_mapb(irot,2)
     do jrot = 1,nrotaa
        q2 = rotaa_mapb(jrot,1)
        q1 = rotaa_mapb(jrot,2)
        AR =  dble(cmat(p2,p1,q1,q2))
        AI = aimag(cmat(p2,p1,q1,q2))
        BR =  dble(cmat(p2,p1,q2,q1))
        BI = aimag(cmat(p2,p1,q2,q1))
        amat(irot,1,jrot,1) = amat(irot,1,jrot,1) + AR - BR
        if (numR == 2) amat(irot,1,jrot,2) = amat(irot,1,jrot,2) - AI - BI
        if (numR == 2) amat(irot,2,jrot,1) = amat(irot,2,jrot,1) + AI - BI
        if (numR == 2) amat(irot,2,jrot,2) = amat(irot,2,jrot,2) + AR + BR
     end do
  end do

  deallocate(cmat)

end subroutine tdcc_xact2_solved_itr_amat_general
!################################################################################
subroutine tdcc_xact2_solved_itr_ainv(nvar, amat, ainv, uvec, vvec)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : iprint
  use mod_control, only : reg_type, throcc3

  implicit none
  integer(c_int), intent(in) :: nvar
  real(kind(0d0)), intent(in) :: amat(1:nvar, 1:nvar)
  real(kind(0d0)), intent(out) :: ainv(1:nvar)
  real(kind(0d0)), intent(out) :: uvec(1:nvar, 1:nvar)
  real(kind(0d0)), intent(out) :: vvec(1:nvar, 1:nvar)
  integer(c_int) :: ivar
  real(kind(0d0)) :: thresh
  real(kind(0d0)) :: inva, denom

  call util_svd_real(nvar, nvar, amat, ainv, uvec, vvec)
  if (iprint > 4) then
     write(6, "('# tdcc_xact2_solved_itr_ainv: singular values of A')")
     do ivar = 1, nvar
        write(6, "(i10, E15.5)") ivar, ainv(ivar)
     end do
  end if

  thresh = throcc3
  do ivar = 1, nvar
     if (reg_type == 0) then      
        ! 1/a <-- a/(a*a + d)
        denom = ainv(ivar) * ainv(ivar) + thresh
        ainv(ivar) = ainv(ivar) / denom
     else if (reg_type == 2) then 
        ! 1/a <-- 1/(a + d * exp(-a / d))
        denom = ainv(ivar) + thresh * exp(-ainv(ivar) / thresh)
        ainv(ivar) = 1d0 / denom
     else
        ! 1/a <-- 1/a
        ainv(ivar) = 1d0 / ainv(ivar)
     end if
  end do

end subroutine tdcc_xact2_solved_itr_ainv
!################################################################################
subroutine tdcc_xact2_solved_itr_bvec(numR, bmath, bmatx, bvec)

  use, intrinsic :: iso_c_binding
  use mod_control, only : icomp
  use mod_ormas, only : nact, nrotaa, rotaa_mapb

  implicit none
  integer(c_int), intent(in) :: numR
  complex(kind(0d0)), intent(in) :: bmath(1:nact, 1:nact)
  complex(kind(0d0)), intent(in) :: bmatx(1:nact, 1:nact)
  real(kind(0d0)), intent(out) :: bvec(1:nrotaa, 1:numR)
  !--------------------------------------------------------------------
  integer(c_int) :: p2, p1, irot
  complex(kind(0d0)) :: Y, Z

  if (icomp == 0) then
     do irot = 1,nrotaa
        p2 = rotaa_mapb(irot,1)
        p1 = rotaa_mapb(irot,2)
        Y = bmath(p2,p1)-conjg(bmath(p1,p2))
        Z = bmatx(p2,p1)-conjg(bmatx(p1,p2))
        bvec(irot,1) = -dble(Y) + dble(Z)
     end do
  else
     do irot = 1,nrotaa
        p2 = rotaa_mapb(irot,1)
        p1 = rotaa_mapb(irot,2)
        Y = bmath(p2,p1)-conjg(bmath(p1,p2))
        Z = bmatx(p2,p1)+conjg(bmatx(p1,p2))
        bvec(irot,1) = aimag(Y) +  dble(Z)
        bvec(irot,2) = -dble(Y) + aimag(Z)
     end do
  end if

end subroutine tdcc_xact2_solved_itr_bvec
!################################################################################
subroutine tdcc_xact2_solved_itr_solve(nvar, ainv, uvec, vvec, bvec, rvec)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : iprint

  implicit none
  integer(c_int), intent(in) :: nvar
  real(kind(0d0)), intent(in) :: ainv(1:nvar)
  real(kind(0d0)), intent(in) :: uvec(1:nvar, 1:nvar)
  real(kind(0d0)), intent(in) :: vvec(1:nvar, 1:nvar)
  real(kind(0d0)), intent(in) :: bvec(1:nvar)
  real(kind(0d0)), intent(out) :: rvec(1:nvar)
  real(kind(0d0)), allocatable :: btmp(:)
  real(kind(0d0)), allocatable :: rtmp(:)
  integer(c_int) :: ivar, jvar

  allocate(btmp(1:nvar))
  allocate(rtmp(1:nvar))

  btmp(1:nvar) = 0d0
  rtmp(1:nvar) = 0d0
  do ivar = 1, nvar
     do jvar = 1, nvar
        btmp(ivar) = btmp(ivar) + uvec(jvar,ivar)*bvec(jvar)
     end do
     rtmp(ivar) = ainv(ivar)*btmp(ivar)
  end do
  if (iprint > 4) then
     do ivar = 1, nvar
        write(6, "('b and U^t*b: ', i5, 4f20.10)") ivar, bvec(ivar), rvec(ivar)
     end do
  end if

  rvec(1:nvar) = 0d0
  do ivar = 1, nvar
     do jvar = 1, nvar
        rvec(ivar) = rvec(ivar) + vvec(ivar,jvar)*rtmp(jvar)
     end do
  end do

  deallocate(rtmp)
  deallocate(btmp)

end subroutine tdcc_xact2_solved_itr_solve
!################################################################################
