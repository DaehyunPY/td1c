!################################################################################
subroutine tdcc_xact2_solveb(bmat, cic, dcic)

  use, intrinsic :: iso_c_binding
  use mod_control, only : icomp
  use mod_const, only : zero,czero,runit,iunit,chalf
  use mod_ormas, only : nact,nrotaa,rotaa_mapb,iprint
  use mod_cc, only : norb1,den1s,den2s,fock,int2x,optbcc
  use mod_cc, only : t1inp,t2inp,g1out

  implicit none
  complex(kind(0d0)), intent(in) :: bmat(1:nact,1:nact)
  complex(kind(0d0)), intent(in) :: cic(1:*)
  complex(kind(0d0)), intent(inout) :: dcic(1:*)
!td1d  complex(kind(0d0)), intent(in) :: cc1(1:norb1,(norb1+1):nact,1:*)
!td1d  complex(kind(0d0)), intent(in) :: cc2((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:*)
!td1d  complex(kind(0d0)), intent(out) :: gcc1(1:norb1,(norb1+1):nact,1:*)
  complex(kind(0d0)) :: tfac,tmp1,tmp2
  integer(c_int) :: iact,jact,aact,bact,pact,qact,ract
  integer(c_int) :: h1,h2,p1,p2,irot,jrot,nrota2
  real(kind(0d0)), allocatable :: l1rhs(:,:)
  real(kind(0d0)), allocatable :: l1vec(:,:)
  real(kind(0d0)), allocatable :: amat(:,:,:,:)

  if (.not.optbcc) stop 'tdcc_xact2_solveb: only optbcc.'

  call tdcc_gettcc1(cic,t1inp)
  call tdcc_gettcc2(cic,t2inp)

  nrota2 = nrotaa*2
  allocate(l1rhs(1:nrotaa,1:2))
  allocate(l1vec(1:nrotaa,1:2))
  allocate(amat(1:nrotaa,1:2,1:nrotaa,1:2))
  l1rhs = zero
  l1vec = zero
  amat = zero
  do irot = 1, nrotaa
     p1 = rotaa_mapb(irot,1)
     h1 = rotaa_mapb(irot,2)
     if (icomp == 1) then
        l1rhs(irot,1) = +aimag(bmat(h1,p1))
        l1rhs(irot,2) = +dble (bmat(h1,p1))
     else
        !GS is assumed to be REAL
        l1rhs(irot,1) = -dble(bmat(h1,p1))
        !GS is assumed to be REAL
     end if
     amat(irot,1,irot,1) = amat(irot,1,irot,1) + runit
     amat(irot,2,irot,2) = amat(irot,2,irot,2) + runit
     do jrot = 1, nrotaa
        p2 = rotaa_mapb(jrot,1) ! particle
        h2 = rotaa_mapb(jrot,2) ! hole

        !##### T1T1 part #####
        if (icomp == 1) then !GS is assumed to be REAL
           amat(irot,1,jrot,2) = amat(irot,1,jrot,2) - aimag(-t1inp(p2,h1,1)*t1inp(p1,h2,1))
           amat(jrot,2,irot,1) = amat(jrot,2,irot,1) - aimag(-t1inp(p2,h1,1)*t1inp(p1,h2,1))
        end if
        if (icomp == 1) then
           amat(irot,1,jrot,1) = amat(irot,1,jrot,1) + dble(-t1inp(p2,h1,1)*t1inp(p1,h2,1))
           amat(irot,2,jrot,2) = amat(irot,2,jrot,2) - dble(-t1inp(p2,h1,1)*t1inp(p1,h2,1))
        else
           amat(irot,1,jrot,1) = amat(irot,1,jrot,1) - dble(-t1inp(p2,h1,1)*t1inp(p1,h2,1))
           amat(irot,2,jrot,2) = amat(irot,2,jrot,2) + dble(-t1inp(p2,h1,1)*t1inp(p1,h2,1))
        end if

        !##### abab part #####
!td1d        if (icomp == 1) then !GS is assumed to be REAL
!td1d           amat(irot,1,jrot,2) = amat(irot,1,jrot,2) - aimag(cc2(p1,p2,h1,h2,2))
!td1d           amat(jrot,2,irot,1) = amat(jrot,2,irot,1) - aimag(cc2(p1,p2,h1,h2,2))
!td1d        end if
!td1d        if (icomp == 1) then
!td1d           amat(irot,1,jrot,1) = amat(irot,1,jrot,1) + dble(cc2(p1,p2,h1,h2,2))
!td1d           amat(irot,2,jrot,2) = amat(irot,2,jrot,2) - dble(cc2(p1,p2,h1,h2,2))
!td1d        else
!td1d           amat(irot,1,jrot,1) = amat(irot,1,jrot,1) - dble(cc2(p1,p2,h1,h2,2))
!td1d           amat(irot,2,jrot,2) = amat(irot,2,jrot,2) + dble(cc2(p1,p2,h1,h2,2))
!td1d        end if

        if (icomp == 1) then !GS is assumed to be REAL
           amat(irot,1,jrot,2) = amat(irot,1,jrot,2) - aimag(t2inp(p1,p2,h1,h2,3))
           amat(jrot,2,irot,1) = amat(jrot,2,irot,1) - aimag(t2inp(p1,p2,h1,h2,3))
        end if
        if (icomp == 1) then
           amat(irot,1,jrot,1) = amat(irot,1,jrot,1) + dble(t2inp(p1,p2,h1,h2,3))
           amat(irot,2,jrot,2) = amat(irot,2,jrot,2) - dble(t2inp(p1,p2,h1,h2,3))
        else
           amat(irot,1,jrot,1) = amat(irot,1,jrot,1) - dble(t2inp(p1,p2,h1,h2,3))
           amat(irot,2,jrot,2) = amat(irot,2,jrot,2) + dble(t2inp(p1,p2,h1,h2,3))
        end if

        !##### aaaa part #####
        if (.not.(p1==p2.or.h1==h2)) then
           if (icomp == 1) then !GS is assumed to be REAL
              if (p1>p2 .and. h1>h2) then
                 amat(irot,1,jrot,2) = amat(irot,1,jrot,2) - aimag(t2inp(p1,p2,h1,h2,1))
                 amat(jrot,2,irot,1) = amat(jrot,2,irot,1) - aimag(t2inp(p1,p2,h1,h2,1))
              else if (p1>p2 .and. h1<h2) then
                 amat(irot,1,jrot,2) = amat(irot,1,jrot,2) + aimag(t2inp(p1,p2,h2,h1,1))
                 amat(jrot,2,irot,1) = amat(jrot,2,irot,1) + aimag(t2inp(p1,p2,h2,h1,1))
              else if (p1<p2 .and. h1>h2) then
                 amat(irot,1,jrot,2) = amat(irot,1,jrot,2) + aimag(t2inp(p2,p1,h1,h2,1))
                 amat(jrot,2,irot,1) = amat(jrot,2,irot,1) + aimag(t2inp(p2,p1,h1,h2,1))
              else if (p1<p2 .and. h1<h2) then
                 amat(irot,1,jrot,2) = amat(irot,1,jrot,2) - aimag(t2inp(p2,p1,h2,h1,1))
                 amat(jrot,2,irot,1) = amat(jrot,2,irot,1) - aimag(t2inp(p2,p1,h2,h1,1))
              end if
           end if
           if (icomp == 1) then
              if (p1>p2 .and. h1>h2) then
                 amat(irot,1,jrot,1) = amat(irot,1,jrot,1) + dble(t2inp(p1,p2,h1,h2,1))
                 amat(irot,2,jrot,2) = amat(irot,2,jrot,2) - dble(t2inp(p1,p2,h1,h2,1))
              else if (p1>p2 .and. h1<h2) then
                 amat(irot,1,jrot,1) = amat(irot,1,jrot,1) - dble(t2inp(p1,p2,h2,h1,1))
                 amat(irot,2,jrot,2) = amat(irot,2,jrot,2) + dble(t2inp(p1,p2,h2,h1,1))
              else if (p1<p2 .and. h1>h2) then
                 amat(irot,1,jrot,1) = amat(irot,1,jrot,1) - dble(t2inp(p2,p1,h1,h2,1))
                 amat(irot,2,jrot,2) = amat(irot,2,jrot,2) + dble(t2inp(p2,p1,h1,h2,1))
              else if (p1<p2 .and. h1<h2) then
                 amat(irot,1,jrot,1) = amat(irot,1,jrot,1) + dble(t2inp(p2,p1,h2,h1,1))
                 amat(irot,2,jrot,2) = amat(irot,2,jrot,2) - dble(t2inp(p2,p1,h2,h1,1))
              end if
           else
              if (p1>p2 .and. h1>h2) then
                 amat(irot,1,jrot,1) = amat(irot,1,jrot,1) - dble(t2inp(p1,p2,h1,h2,1))
                 amat(irot,2,jrot,2) = amat(irot,2,jrot,2) + dble(t2inp(p1,p2,h1,h2,1))
              else if (p1>p2 .and. h1<h2) then
                 amat(irot,1,jrot,1) = amat(irot,1,jrot,1) + dble(t2inp(p1,p2,h2,h1,1))
                 amat(irot,2,jrot,2) = amat(irot,2,jrot,2) - dble(t2inp(p1,p2,h2,h1,1))
              else if (p1<p2 .and. h1>h2) then
                 amat(irot,1,jrot,1) = amat(irot,1,jrot,1) + dble(t2inp(p2,p1,h1,h2,1))
                 amat(irot,2,jrot,2) = amat(irot,2,jrot,2) - dble(t2inp(p2,p1,h1,h2,1))
              else if (p1<p2 .and. h1<h2) then
                 amat(irot,1,jrot,1) = amat(irot,1,jrot,1) - dble(t2inp(p2,p1,h2,h1,1))
                 amat(irot,2,jrot,2) = amat(irot,2,jrot,2) + dble(t2inp(p2,p1,h2,h1,1))
              end if
           end if
        end if
     end do
  end do

!     if (iprint > 4) then
!        write(6, "('# ccnew_rmat2_bccd_l1vec: amat-11 ')")
!        do irot = 1, nrotaa
!           do jrot = 1, nrotaa
!              write(6, "(e12.4)", advance='no') amat(irot,1,jrot,1)
!           end do
!           write(6, *)
!        end do
!        write(6, "('# ccnew_rmat2_bccd_l1vec: amat-12 ')")
!        do irot = 1, nrotaa
!           do jrot = 1, nrotaa
!              write(6, "(e12.4)", advance='no') amat(irot,1,jrot,2)
!           end do
!           write(6, *)
!        end do
!        write(6, "('# ccnew_rmat2_bccd_l1vec: amat-21 ')")
!        do irot = 1, nrotaa
!           do jrot = 1, nrotaa
!              write(6, "(e12.4)", advance='no') amat(irot,2,jrot,1)
!           end do
!           write(6, *)
!        end do
!        write(6, "('# ccnew_rmat2_bccd_l1vec: amat-22 ')")
!        do irot = 1, nrotaa
!           do jrot = 1, nrotaa
!              write(6, "(e12.4)", advance='no') amat(irot,2,jrot,2)
!           end do
!           write(6, *)
!        end do
!     end if

  call tdcc_bcct1_lineq(nrota2, amat, l1rhs, l1vec)
! call ccnew_rmat2_bccd_lineq(-1, nrota2, amat, l1rhs, l1vec)
! call ccnew_rmat2_bccd_lineq2(abs(dpsi_reg), nrota2, amat, l1rhs, l1vec)
! call ccnew_rmat2_bccd_lineq2(-1, nrota2, amat, l1rhs, l1vec)

  g1out = 0d0
  do irot = 1, nrotaa
     p1 = rotaa_mapb(irot,1)
     h1 = rotaa_mapb(irot,2)
     if (icomp == 1) then
        g1out(h1,p1,1) = l1vec(irot,1)+iunit*l1vec(irot,2)
     else
        !GS is assumed to be REAL
        g1out(h1,p1,1) = l1vec(irot,1)
        !GS is assumed to be REAL
     end if
  end do

! <<<<< HACK FOR TDCC_SCALE
  if (icomp == 1) then
     tfac = -iunit
  else
     tfac = -runit
  end if
  g1out(:,:,1) = g1out(:,:,1)*tfac
  g1out(:,:,2) = g1out(:,:,1)
! <<<<< HACK FOR TDCC_SCALE

  call tdcc_putgcc1(g1out,dcic)

  deallocate(amat)
  deallocate(l1vec)
  deallocate(l1rhs)

end subroutine tdcc_xact2_solveb
!################################################################################
