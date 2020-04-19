!######################################################################
subroutine tdcc_bcct1(cic,xmat)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only :nact
  use mod_cc, only : norb1,cc_rank,bcc,optbcc,t2inp,t3inp,t1out

  implicit none
  complex(kind(0d0)), intent(in) :: cic(1:*)
  complex(kind(0d0)), intent(out) :: xmat(1:nact,1:nact)

  if (.not.(bcc.or.optbcc)) return
  t1out = 0d0

  if (cc_rank == 2) then
     call tdcc_gettcc2(cic,t2inp)
     call bccd_t1p_main
    !call ccsd_t1_main
  end if
  if (cc_rank == 3) then
     call tdcc_gettcc2(cic,t2inp)
     call tdcc_gettcc3(cic,t3inp)
     call bccdt_t1p_main
  end if

  call tdcc_gettcc2(cic,t2inp)
  call tdcc_bcct1_xmat(t1out,t2inp,xmat)

end subroutine tdcc_bcct1
!################################################################################
subroutine tdcc_bcct1_xmat(hcc1, cc2, xmat)

  use, intrinsic :: iso_c_binding
  use mod_control, only : icomp
  use mod_const, only : zero,runit,iunit
  use mod_ormas, only : nact,nrotaa,rotaa_mapb
  use mod_cc, only : norb1

  implicit none
  complex(kind(0d0)), intent(in) :: hcc1((norb1+1):nact,1:norb1,1:*)
  complex(kind(0d0)), intent(in) :: cc2((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:6)
  complex(kind(0d0)), intent(out) :: xmat(1:nact,1:nact)
  integer(c_int) :: h1,h2,p1,p2,irot,jrot,nrota2
  real(kind(0d0)), allocatable :: x1rhs(:,:)
  real(kind(0d0)), allocatable :: x1vec(:,:)
  real(kind(0d0)), allocatable :: amat(:,:,:,:)

!DEBUG
!  do p1 = norb1+1, nact
!  do h1 = 1, norb1
!     xmat(p1,h1) = hcc1(p1,h1,1)
!     xmat(h1,p1) = -conjg(xmat(p1,h1))
!  end do
!  end do
!  return
!DEBUG

  nrota2 = nrotaa*2
  allocate(x1rhs(1:nrotaa,1:2))
  allocate(x1vec(1:nrotaa,1:2))
  allocate(amat(1:nrotaa,1:2,1:nrotaa,1:2))
  x1rhs = zero
  x1vec = zero
  amat = zero
  do irot = 1, nrotaa
     p1 = rotaa_mapb(irot,1)
     h1 = rotaa_mapb(irot,2)
     if (icomp == 1) then
        x1rhs(irot,1) = aimag(hcc1(p1,h1,1))
        x1rhs(irot,2) = -dble(hcc1(p1,h1,1))
     else
        !GS is assumed to be REAL
        x1rhs(irot,1) = -dble (hcc1(p1,h1,1))
        !GS is assumed to be REAL
     end if
     amat(irot,1,irot,1) = amat(irot,1,irot,1) + runit
     amat(irot,2,irot,2) = amat(irot,2,irot,2) + runit
     do jrot = 1, nrotaa
        p2 = rotaa_mapb(jrot,1) ! particle
        h2 = rotaa_mapb(jrot,2) ! hole

        !##### abab part #####        
!td1d        if (icomp == 1) then !GS is assumed to be REAL
!td1d           amat(irot,1,jrot,2) = amat(irot,1,jrot,2) - aimag(cc2(p1,p2,h1,h2,2))
!td1d           amat(jrot,2,irot,1) = amat(jrot,2,irot,1) - aimag(cc2(p1,p2,h1,h2,2))
!td1d        end if
!td1d        amat(irot,1,jrot,1) = amat(irot,1,jrot,1) - dble(cc2(p1,p2,h1,h2,2))
!td1d        amat(irot,2,jrot,2) = amat(irot,2,jrot,2) + dble(cc2(p1,p2,h1,h2,2))
        if (icomp == 1) then !GS is assumed to be REAL
           amat(irot,1,jrot,2) = amat(irot,1,jrot,2) - aimag(cc2(p1,p2,h1,h2,3))
           amat(jrot,2,irot,1) = amat(jrot,2,irot,1) - aimag(cc2(p1,p2,h1,h2,3))
        end if
        amat(irot,1,jrot,1) = amat(irot,1,jrot,1) - dble(cc2(p1,p2,h1,h2,3))
        amat(irot,2,jrot,2) = amat(irot,2,jrot,2) + dble(cc2(p1,p2,h1,h2,3))

        !##### aaaa part #####
        if (.not.(p1==p2.or.h1==h2)) then
           if (icomp == 1) then !GS is assumed to be REAL
              if (p1>p2 .and. h1>h2) then
                 amat(irot,1,jrot,2) = amat(irot,1,jrot,2) - aimag(cc2(p1,p2,h1,h2,1))
                 amat(jrot,2,irot,1) = amat(jrot,2,irot,1) - aimag(cc2(p1,p2,h1,h2,1))
              else if (p1>p2 .and. h1<h2) then
                 amat(irot,1,jrot,2) = amat(irot,1,jrot,2) + aimag(cc2(p1,p2,h2,h1,1))
                 amat(jrot,2,irot,1) = amat(jrot,2,irot,1) + aimag(cc2(p1,p2,h2,h1,1))
              else if (p1<p2 .and. h1>h2) then
                 amat(irot,1,jrot,2) = amat(irot,1,jrot,2) + aimag(cc2(p2,p1,h1,h2,1))
                 amat(jrot,2,irot,1) = amat(jrot,2,irot,1) + aimag(cc2(p2,p1,h1,h2,1))
              else if (p1<p2 .and. h1<h2) then
                 amat(irot,1,jrot,2) = amat(irot,1,jrot,2) - aimag(cc2(p2,p1,h2,h1,1))
                 amat(jrot,2,irot,1) = amat(jrot,2,irot,1) - aimag(cc2(p2,p1,h2,h1,1))
              end if
           end if
           if (p1>p2 .and. h1>h2) then
              amat(irot,1,jrot,1) = amat(irot,1,jrot,1) - dble(cc2(p1,p2,h1,h2,1))
              amat(irot,2,jrot,2) = amat(irot,2,jrot,2) + dble(cc2(p1,p2,h1,h2,1))
           else if (p1>p2 .and. h1<h2) then
              amat(irot,1,jrot,1) = amat(irot,1,jrot,1) + dble(cc2(p1,p2,h2,h1,1))
              amat(irot,2,jrot,2) = amat(irot,2,jrot,2) - dble(cc2(p1,p2,h2,h1,1))
           else if (p1<p2 .and. h1>h2) then
              amat(irot,1,jrot,1) = amat(irot,1,jrot,1) + dble(cc2(p2,p1,h1,h2,1))
              amat(irot,2,jrot,2) = amat(irot,2,jrot,2) - dble(cc2(p2,p1,h1,h2,1))
           else if (p1<p2 .and. h1<h2) then
              amat(irot,1,jrot,1) = amat(irot,1,jrot,1) - dble(cc2(p2,p1,h2,h1,1))
              amat(irot,2,jrot,2) = amat(irot,2,jrot,2) + dble(cc2(p2,p1,h2,h1,1))
           end if
        end if
     end do
  end do

!     if (iprint > 4) then
!        write(6, "('# tdcc_rmat2_bccd_x1vec: amat-11 ')")
!        do irot = 1, nrotaa
!           do jrot = 1, nrotaa
!              write(6, "(e12.4)", advance='no') amat(irot,1,jrot,1)
!           end do
!           write(6, *)
!        end do
!        write(6, "('# tdcc_rmat2_bccd_x1vec: amat-12 ')")
!        do irot = 1, nrotaa
!           do jrot = 1, nrotaa
!              write(6, "(e12.4)", advance='no') amat(irot,1,jrot,2)
!           end do
!           write(6, *)
!        end do
!        write(6, "('# tdcc_rmat2_bccd_x1vec: amat-21 ')")
!        do irot = 1, nrotaa
!           do jrot = 1, nrotaa
!              write(6, "(e12.4)", advance='no') amat(irot,2,jrot,1)
!           end do
!           write(6, *)
!        end do
!        write(6, "('# tdcc_rmat2_bccd_x1vec: amat-22 ')")
!        do irot = 1, nrotaa
!           do jrot = 1, nrotaa
!              write(6, "(e12.4)", advance='no') amat(irot,2,jrot,2)
!           end do
!           write(6, *)
!        end do
!     end if

!OLD  call tdcc_rmat2_bccd_lineq(-1, nrota2, amat, x1rhs, x1vec)
  call tdcc_bcct1_lineq(nrota2, amat, x1rhs, x1vec)

  do irot = 1, nrotaa
     p1 = rotaa_mapb(irot,1)
     h1 = rotaa_mapb(irot,2)
     if (icomp == 1) then
        xmat(p1,h1) = x1vec(irot,1)+iunit*x1vec(irot,2)
        xmat(h1,p1) = -conjg(xmat(p1,h1))
     else
        !GS is assumed to be REAL
        xmat(p1,h1) = x1vec(irot,1)
        xmat(h1,p1) = -x1vec(irot,1)
        !GS is assumed to be REAL
     end if
  end do

  deallocate(amat)
  deallocate(x1vec)
  deallocate(x1rhs)

end subroutine tdcc_bcct1_xmat
!################################################################################
subroutine tdcc_bcct1_lineq(nvar, amat, bvec, rvec)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : iprint
  use mod_const, only : one, czero
  use mod_control, only : reg_type, throcc3

  implicit none
  integer(c_int), intent(in) :: nvar
  real(kind(0d0)), intent(in) :: amat(1:nvar, 1:nvar)
  real(kind(0d0)), intent(in) :: bvec(1:nvar)
  real(kind(0d0)), intent(out) :: rvec(1:nvar)

  integer(c_int) :: ivar, jvar
  real(kind(0d0)) :: thresh, inva, denom
  real(kind(0d0)), allocatable :: aeig(:)
  real(kind(0d0)), allocatable :: uvec(:,:), vvec(:,:), btmp(:), rtmp(:)

  thresh = throcc3
  allocate(aeig(1:nvar))
  allocate(uvec(1:nvar, 1:nvar))
  allocate(vvec(1:nvar, 1:nvar))
  allocate(btmp(1:nvar))
  allocate(rtmp(1:nvar))

!OLD  call util_svd_real(nvar, nvar, amat, aeig, uvec, vvec)
  call util_diag_real(.true., nvar, amat, uvec)
  vvec = uvec
  do ivar = 1, nvar
     aeig(ivar) = amat(ivar,ivar)
  end do
  !dpsi_smin = aeig(nvar)

  if (iprint > 4) then
     write(6, "('# tdcc_rmat2_lineq: singular values of A')")
     do ivar = 1, min(8, nvar)
        write(6, "(E15.5)", advance = 'no') aeig(ivar)
     end do
     write(6,*)
  end if

  btmp(1:nvar) = czero
  do ivar = 1, nvar
     do jvar = 1, nvar
        btmp(ivar) = btmp(ivar) + uvec(jvar, ivar) * bvec(jvar)
     end do
  end do

!  if (iprint > 1) then
!     do ivar = 1, nvar
!        write(6, "('b and U^t*b: ', i5, 4f20.10)") ivar, bvec(ivar), btmp(ivar)
!     end do
!  end if

  rtmp(1:nvar) = czero
  do ivar = 1, nvar
     if (reg_type == 0) then      ! 1/a <-- a/(a*a + d*d)
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
        stop 'bad reg_type in tdcc_dpsidt_donly_solv.'
     end if
     rtmp(ivar) = inva * btmp(ivar)
  end do

  rvec(1:nvar) = czero
  do ivar = 1, nvar
     do jvar = 1, nvar
        rvec(ivar) = rvec(ivar) + vvec(ivar, jvar) * rtmp(jvar)
     end do
  end do

  deallocate(rtmp)
  deallocate(btmp)
  deallocate(vvec)
  deallocate(uvec)
  deallocate(aeig)

end subroutine tdcc_bcct1_lineq
!######################################################################
