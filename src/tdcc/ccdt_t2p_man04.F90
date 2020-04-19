!##########################################################
subroutine ccdt_t2p_man04(i0,work1,work2,work3)

!1: i0 ( a b i j )_v + = 1 * v ( a b i j )_v 0
!6: i0 ( a b i j )_vt + = 1/2 * Sum ( c d ) * t ( c d i j )_t * v ( a b c d )_v 0
!7: i0 ( a b i j )_tf + = 1 * Sum ( c k ) * t ( a b c i j k )_t * f ( k c )_f 0
!8: i0 ( a b i j )_vt + = -1/2 * P( 2 ) * Sum ( k l c ) * t ( a b c i k l )_t * v ( k l j c )_v 0
!9: i0 ( a b i j )_vt + = -1/2 * P( 2 ) * Sum ( k c d ) * t ( a c d i j k )_t * v ( k b c d )_v 0

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact,nact1,act1_ll,act1_ul
  use mod_cc,only : fock,int2x,norb1,ncc2aa,ncc2ab,t2inp,t3inp,cc_rank
  use mod_cc,only : h1_cc2aa,h2_cc2aa,p1_cc2aa,p2_cc2aa
  use mod_cc,only : h1_cc2ab,h2_cc2ab,p1_cc2ab,p2_cc2ab
  use mod_cc,only : ncc2aa_act1,ncc2ab_act1
  use mod_cc,only : h1_cc2aa_act1,h2_cc2aa_act1,p1_cc2aa_act1,p2_cc2aa_act1
  use mod_cc,only : h1_cc2ab_act1,h2_cc2ab_act1,p1_cc2ab_act1,p2_cc2ab_act1
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i0((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:2)
  complex(kind(0d0)),intent(inout) :: work1(1),work2(1),work3(1)
  complex(kind(0d0)) :: tmp,tmp1,tmp2

  integer(c_int) :: icc,a,b,c,d,e,i,j,k,l,m

  !$omp parallel default(shared) private(a,b,i,j,tmp,tmp1,tmp2)
  ! ##### CCD contributions #####
  !$omp do
  do icc = 1,ncc2aa
     a = p1_cc2aa(icc)
     b = p2_cc2aa(icc)
     i = h1_cc2aa(icc)
     j = h2_cc2aa(icc)
     tmp = 0d0
     ! diagram 1
     tmp = tmp + int2x(a,b,i,j,spin_int2aa)
     ! diagram 6
     do c = norb1+1,nact
     do d = norb1+1,c-1
        tmp = tmp + t2inp(c,d,i,j,spin_t2aa)*int2x(a,b,c,d,spin_int2aa)
     end do
     end do
     i0(a,b,i,j,1) = i0(a,b,i,j,1) + tmp
  end do
  !$omp end do
  !$omp do
  do icc = 1,ncc2ab
     a = p1_cc2ab(icc)
     b = p2_cc2ab(icc)
     i = h1_cc2ab(icc)
     j = h2_cc2ab(icc)
     tmp = 0d0
     ! diagram 1
     tmp = tmp + int2x(a,b,i,j,spin_int2ab)
     ! diagram 6
     do c = norb1+1,nact
     do d = norb1+1,nact
        tmp = tmp + t2inp(c,d,i,j,spin_t2ab)*int2x(a,b,c,d,spin_int2ab)
     end do
     end do
     i0(a,b,i,j,2) = i0(a,b,i,j,2) + tmp
  end do
  !$omp end do

  ! ##### CCDT contributions #####
  if (cc_rank >= 3) then
     !$omp do
     do icc = 1,ncc2aa_act1
        a = p1_cc2aa_act1(icc)
        b = p2_cc2aa_act1(icc)
        i = h1_cc2aa_act1(icc)
        j = h2_cc2aa_act1(icc)
        tmp = 0d0
        ! diagram 7
        do c = norb1+1,act1_ul
        do k = act1_ll,norb1
           tmp = tmp + (t3inp(a,b,c,i,j,k,spin_t3aaa) &
                       +t3inp(a,b,c,i,j,k,spin_t3aab))*fock(k,c,spin_focka)
        end do
        end do
        ! diagram 8
        do c = norb1+1,act1_ul
        do k = act1_ll,norb1
           do l = act1_ll,k-1
              tmp = tmp - t3inp(a,b,c,i,k,l,spin_t3aaa)*int2x(k,l,j,c,spin_int2aa) &
                        + t3inp(a,b,c,j,k,l,spin_t3aaa)*int2x(k,l,i,c,spin_int2aa)
           end do
           do l = act1_ll,norb1
              tmp = tmp - t3inp(a,b,c,i,k,l,spin_t3aab)*int2x(k,l,j,c,spin_int2ab) &
                        + t3inp(a,b,c,j,k,l,spin_t3aab)*int2x(k,l,i,c,spin_int2ab)
           end do
        end do
        end do
        ! diagram 9
        do k = act1_ll,norb1
        do c = norb1+1,act1_ul
           do d = norb1+1,c-1
              tmp = tmp - t3inp(a,c,d,i,j,k,spin_t3aaa)*int2x(k,b,c,d,spin_int2aa) &
                        + t3inp(b,c,d,i,j,k,spin_t3aaa)*int2x(k,a,c,d,spin_int2aa)
           end do
           do d = norb1+1,act1_ul
              tmp = tmp + t3inp(a,c,d,i,j,k,spin_t3aab)*int2x(b,k,c,d,spin_int2ab) &
                        - t3inp(b,c,d,i,j,k,spin_t3aab)*int2x(a,k,c,d,spin_int2ab)
           end do
        end do
        end do
        i0(a,b,i,j,1) = i0(a,b,i,j,1) + tmp
     end do
     !$omp end do
     !$omp do
     do icc = 1,ncc2ab_act1
        a = p1_cc2ab_act1(icc)
        b = p2_cc2ab_act1(icc)
        i = h1_cc2ab_act1(icc)
        j = h2_cc2ab_act1(icc)
        tmp = 0d0
        ! diagram 7
        do c = norb1+1,act1_ul
        do k = act1_ll,norb1
           tmp = tmp + (t3inp(a,c,b,i,k,j,spin_t3aab) &
                       +t3inp(b,c,a,j,k,i,spin_t3aab))*fock(k,c,spin_focka)
        end do
        end do
        ! diagram 8
        do c = norb1+1,act1_ul
        do k = act1_ll,norb1
           do l = act1_ll,k-1
              tmp = tmp - t3inp(b,c,a,k,l,i,spin_t3aab)*int2x(k,l,j,c,spin_int2aa) &
                        + t3inp(a,c,b,l,k,j,spin_t3aab)*int2x(k,l,i,c,spin_int2aa)
           end do
           do l = act1_ll,norb1
              tmp = tmp - t3inp(a,c,b,i,l,k,spin_t3aab)*int2x(k,l,j,c,spin_int2ab) &
                        + t3inp(c,b,a,j,l,k,spin_t3aab)*int2x(k,l,i,c,spin_int2ab)
           end do
        end do
        end do
        ! diagram 9
        do k = act1_ll,norb1
        do c = norb1+1,act1_ul
           do d = norb1+1,c-1
              tmp = tmp - t3inp(c,d,a,j,k,i,spin_t3aab)*int2x(k,b,c,d,spin_int2aa) &
                        + t3inp(d,c,b,i,k,j,spin_t3aab)*int2x(k,a,c,d,spin_int2aa)
           end do
           do d = norb1+1,act1_ul
              tmp = tmp + t3inp(a,c,d,i,k,j,spin_t3aab)*int2x(k,b,c,d,spin_int2ab) &
                        - t3inp(b,d,c,k,j,i,spin_t3aab)*int2x(a,k,c,d,spin_int2ab)
           end do
        end do
        end do
        i0(a,b,i,j,2) = i0(a,b,i,j,2) + tmp
     end do
     !$omp end do

     !if (.false.) then
     if (nact1 < nact) then
        ! diagram 8aa
        !$omp do collapse(3)
        do a = norb1+1,act1_ul
        do i = 1,act1_ll-1
        do j = act1_ll,norb1
        do b = a+1,act1_ul
           tmp = 0d0
           do c = norb1+1,act1_ul
           do k = act1_ll,norb1
              do l = act1_ll,k-1
                 tmp = tmp + t3inp(a,b,c,j,k,l,spin_t3aaa)*int2x(k,l,i,c,spin_int2aa)
              end do
              do l = act1_ll,norb1
                 tmp = tmp + t3inp(a,b,c,j,k,l,spin_t3aab)*int2x(k,l,i,c,spin_int2ab)
              end do
           end do
           end do
           i0(a,b,i,j,1) = i0(a,b,i,j,1) + tmp
        end do
        end do
        end do
        end do
        !$omp end do

        ! diagram 9aa
        !$omp do collapse(3)
        do a = act1_ul+1,nact
        do b = norb1+1,act1_ul
        do i = act1_ll,norb1
        do j = i+1,norb1
           do k = act1_ll,norb1
           do c = norb1+1,act1_ul
              do d = norb1+1,c-1
                 tmp = tmp + t3inp(b,c,d,i,j,k,spin_t3aaa)*int2x(k,a,c,d,spin_int2aa)
              end do
              do d = norb1+1,act1_ul
                 tmp = tmp - t3inp(b,c,d,i,j,k,spin_t3aab)*int2x(a,k,c,d,spin_int2ab)
              end do
           end do
           end do
           i0(a,b,i,j,1) = i0(a,b,i,j,1) + tmp
        end do
        end do
        end do
        end do

        ! diagram 8ab
        !$omp do collapse(4)
        do i = 1,act1_ll-1
        do j = act1_ll,norb1
        do a = norb1+1,act1_ul
        do b = norb1+1,act1_ul
           tmp1 = 0d0
           tmp2 = 0d0
           do c = norb1+1,act1_ul
           do k = act1_ll,norb1
              do l = act1_ll,k-1
                 tmp1 = tmp1 - t3inp(b,c,a,k,l,j,spin_t3aab)*int2x(k,l,i,c,spin_int2aa)
                 tmp2 = tmp2 + t3inp(a,c,b,l,k,j,spin_t3aab)*int2x(k,l,i,c,spin_int2aa)
              end do
              do l = act1_ll,norb1
                 tmp1 = tmp1 - t3inp(a,c,b,j,l,k,spin_t3aab)*int2x(k,l,i,c,spin_int2ab)
                 tmp2 = tmp2 + t3inp(c,b,a,j,l,k,spin_t3aab)*int2x(k,l,i,c,spin_int2ab)
              end do
           end do
           end do
           i0(a,b,j,i,2) = i0(a,b,j,i,2) + tmp1
           i0(a,b,i,j,2) = i0(a,b,i,j,2) + tmp2
        end do
        end do
        end do
        end do
        !$omp end do

        ! diagram 9ab
        !$omp do collapse(4)
        do a = act1_ul+1,nact
        do b = norb1+1,act1_ul
        do i = 1,act1_ll-1
        do j = act1_ll,norb1
           tmp1 = 0d0
           tmp2 = 0d0
           do k = act1_ll,norb1
           do c = norb1+1,act1_ul
              do d = norb1+1,c-1
                 tmp1 = tmp1 - t3inp(c,d,b,j,k,i,spin_t3aab)*int2x(k,a,c,d,spin_int2aa)
                 tmp2 = tmp2 + t3inp(d,c,b,i,k,j,spin_t3aab)*int2x(k,a,c,d,spin_int2aa)
              end do
              do d = norb1+1,act1_ul
                 tmp1 = tmp1 + t3inp(b,c,d,i,k,j,spin_t3aab)*int2x(k,a,c,d,spin_int2ab)
                 tmp2 = tmp2 - t3inp(b,d,c,k,j,i,spin_t3aab)*int2x(a,k,c,d,spin_int2ab)
              end do
           end do
           end do
           i0(b,a,i,j,2) = i0(b,a,i,j,2) + tmp1
           i0(a,b,i,j,2) = i0(a,b,i,j,2) + tmp2
        end do
        end do
        end do
        end do
        !$omp end do
     end if
  end if
  !$omp end parallel

end subroutine ccdt_t2p_man04
!##########################################################
