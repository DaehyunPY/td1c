!**********************************************************
subroutine ccdt_bmat2_main(fac, bmat)

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact,rotaa_mapb
  use mod_cc,only : cc_code
  use mod_cc,only : norb1,t2inp,g2inp,dt2inp,dg2inp
  use mod_cc,only : nXai,nXij,nXab,aiX,ijX,abX
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(in) :: fac
  complex(kind(0d0)),intent(inout) :: bmat(1:nact,1:nact)
  complex(kind(0d0)),allocatable :: b2(:,:),b3(:,:)
  integer(c_int) :: a,b,i,j,ai,i1,i2,a1,a2,i1i2,a1a2,irot
  complex(kind(0d0)) :: tmp1,tmp2

  allocate(b2(nact,nact))
  allocate(b3(nact,nact))
  b2 = 0d0
  b3 = 0d0
  if (cc_code(1:4) == 'auto') then
     call ccdt_bmat2p_1(1,1,b2)
     call ccdt_bmat2p_2(1,1,b2)
     call ccdt_bmat2p_3(1,1,b2)
     call ccdt_bmat3p_1(1,1,b3)
     call ccdt_bmat3p_2(1,1,b3)
  else
     call ccdt_bmat2p_man01(b2,cc_work1,cc_work2)
     call ccdt_bmat2p_man02(b2,cc_work1,cc_work2)
     call ccdt_bmat2p_man03(b2,cc_work1,cc_work2)
     call ccdt_bmat3p_man01(b3,cc_work1,cc_work2)
     call ccdt_bmat3p_man02(b3,cc_work1,cc_work2)
  end if
  do ai = 1, nXai
     irot = aiX(ai)
     a = rotaa_mapb(irot,1)
     i = rotaa_mapb(irot,2)
     bmat(a,i) = bmat(a,i) + fac*(b2(i,a) - b3(i,a))
  end do
  deallocate(b2)
  deallocate(b3)

!DEBUG
!return
!DEBUG

  do i1i2 = 1, nXij
     irot = ijX(i1i2)
     i1 = rotaa_mapb(irot,1)
     i2 = rotaa_mapb(irot,2)
     tmp1 = 0d0
     tmp2 = 0d0
     do j = 1, norb1
     do a = norb1+1, nact
        do b = norb1+1, a-1
           tmp1 = tmp1 -  g2inp(i1,j,a,b,spin_g2aa)*dt2inp(a,b,i2,j,spin_t2aa) &
                       + dg2inp(i1,j,a,b,spin_g2aa)* t2inp(a,b,i2,j,spin_t2aa)
           tmp2 = tmp2 -  g2inp(i2,j,a,b,spin_g2aa)*dt2inp(a,b,i1,j,spin_t2aa) &
                       + dg2inp(i2,j,a,b,spin_g2aa)* t2inp(a,b,i1,j,spin_t2aa)
        end do
        do b = norb1+1, nact
           tmp1 = tmp1 -  g2inp(i1,j,a,b,spin_g2ab)*dt2inp(a,b,i2,j,spin_t2ab) &
                       + dg2inp(i1,j,a,b,spin_g2ab)* t2inp(a,b,i2,j,spin_t2ab)
           tmp2 = tmp2 -  g2inp(i2,j,a,b,spin_g2ab)*dt2inp(a,b,i1,j,spin_t2ab) &
                       + dg2inp(i2,j,a,b,spin_g2ab)* t2inp(a,b,i1,j,spin_t2ab)
        end do
     end do
     end do
     bmat(i1,i2) = bmat(i1,i2) + fac * tmp1
     bmat(i2,i1) = bmat(i2,i1) + fac * tmp2
  end do

  do a1a2 = 1, nXab
     irot = abX(a1a2)
     a2 = rotaa_mapb(irot,1)
     a1 = rotaa_mapb(irot,2)
     tmp1 = 0d0
     tmp2 = 0d0
     do b = norb1+1, nact
     do i = 1, norb1
        do j = 1, i-1
           tmp1 = tmp1 +  g2inp(i,j,a1,b,spin_g2aa)*dt2inp(a2,b,i,j,spin_t2aa) &
                       - dg2inp(i,j,a1,b,spin_g2aa)* t2inp(a2,b,i,j,spin_t2aa)
           tmp2 = tmp2 +  g2inp(i,j,a2,b,spin_g2aa)*dt2inp(a1,b,i,j,spin_t2aa) &
                       - dg2inp(i,j,a2,b,spin_g2aa)* t2inp(a1,b,i,j,spin_t2aa)
        end do
        do j = 1, norb1
           tmp1 = tmp1 +  g2inp(i,j,a1,b,spin_g2ab)*dt2inp(a2,b,i,j,spin_t2ab) &
                       - dg2inp(i,j,a1,b,spin_g2ab)* t2inp(a2,b,i,j,spin_t2ab)
           tmp2 = tmp2 +  g2inp(i,j,a2,b,spin_g2ab)*dt2inp(a1,b,i,j,spin_t2ab) &
                       - dg2inp(i,j,a2,b,spin_g2ab)* t2inp(a1,b,i,j,spin_t2ab)
        end do
     end do
     end do
     bmat(a2,a1) = bmat(a2,a1) + fac * tmp1
     bmat(a1,a2) = bmat(a1,a2) + fac * tmp2
  end do

end subroutine ccdt_bmat2_main
!**********************************************************
