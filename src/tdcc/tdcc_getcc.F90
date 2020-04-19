!######################################################################
subroutine tdcc_getcc(cic)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nact
  use mod_cc, only : cc_rank,t0inp,g0inp,t1inp,g1inp,t2inp,g2inp,t3inp,g3inp
  implicit none
  complex(c_double_complex), intent(in) :: cic(1:*)

  call tdcc_gettcc0(cic,t0inp)
  if (cc_rank >= 1) call tdcc_gettcc1(cic,t1inp)
  if (cc_rank >= 2) call tdcc_gettcc2(cic,t2inp)
  if (cc_rank >= 3) call tdcc_gettcc3(cic,t3inp)

  call tdcc_getgcc0(cic,g0inp)
  if (cc_rank >= 1) call tdcc_getgcc1(cic,g1inp)
  if (cc_rank >= 2) call tdcc_getgcc2(cic,g2inp)
  if (cc_rank >= 3) call tdcc_getgcc3(cic,g3inp)

end subroutine tdcc_getcc
!######################################################################
subroutine tdcc_gettcc0(cic, tcc0)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : ndetx
  use mod_cc, only : map_cc0
  implicit none
  complex(c_double_complex), intent(in) :: cic(1:ndetx,1:2)
  complex(c_double_complex), intent(out) :: tcc0

  tcc0 = cic(map_cc0(1,3),1)

end subroutine tdcc_gettcc0
!######################################################################
subroutine tdcc_getgcc0(cic, gcc0)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : ndetx
  use mod_cc, only : map_cc0
  implicit none
  complex(c_double_complex), intent(in) :: cic(1:ndetx,1:2)
  complex(c_double_complex), intent(out) :: gcc0

  gcc0 = cic(map_cc0(1,3),2)

end subroutine tdcc_getgcc0
!######################################################################
subroutine tdcc_gettcc1(cic, tcc1)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nact,ndetx,nelact
  use mod_cc, only : norb1,ncc1a
  use mod_cc, only : h1_cc1a,p1_cc1a,map_cc1a
  implicit none
  complex(c_double_complex), intent(in) :: cic(1:ndetx,1:2)
  complex(c_double_complex), intent(out) :: tcc1((norb1+1):nact,1:norb1,1:2)
  integer(c_int) :: icc,p1,h1,idet

  tcc1 = 0d0
  if (nelact(3)< 1) return
  do icc = 1, ncc1a
     p1 = p1_cc1a(icc)
     h1 = h1_cc1a(icc)
     idet = map_cc1a(icc,3)
     tcc1(p1,h1,1) = cic(idet,1)
  end do
  tcc1(:,:,2) = tcc1(:,:,1)

end subroutine tdcc_gettcc1
!######################################################################
subroutine tdcc_getgcc1(cic, gcc1)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nact,ndetx,nelact
  use mod_cc, only : norb1,ncc1a
  use mod_cc, only : h1_cc1a,p1_cc1a,map_cc1a
  implicit none
  complex(c_double_complex), intent(in) :: cic(1:ndetx,1:2)
  complex(c_double_complex), intent(out) :: gcc1(1:norb1,(norb1+1):nact,1:2)
  integer(c_int) :: icc,p1,h1,idet

  gcc1 = 0d0
  if (nelact(3)< 1) return
  do icc = 1, ncc1a
     p1 = p1_cc1a(icc)
     h1 = h1_cc1a(icc)
     idet = map_cc1a(icc,3)
     gcc1(h1,p1,1) = cic(idet,2)
  end do
  gcc1(:,:,2) = gcc1(:,:,1)

end subroutine tdcc_getgcc1
!######################################################################
subroutine tdcc_gettcc2(cic, tcc2)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nact,ndetx,nelact
  use mod_cc, only : norb1,ncc2aa,ncc2ab,map_cc2aa,map_cc2ab
  use mod_cc, only : h1_cc2aa,h2_cc2aa,p1_cc2aa,p2_cc2aa
  use mod_cc, only : h1_cc2ab,h2_cc2ab,p1_cc2ab,p2_cc2ab
  implicit none
  complex(c_double_complex), intent(in) :: cic(1:ndetx,1:2)
  complex(c_double_complex), intent(out) :: &
       tcc2((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:6)
  integer(c_int) :: icc,p1,p2,h1,h2,idet

  tcc2 = 0d0
  if (nelact(3)< 2) return
  if (nelact(1)>=2) then
     do icc = 1, ncc2aa
        p1 = p1_cc2aa(icc)
        p2 = p2_cc2aa(icc)
        h1 = h1_cc2aa(icc)
        h2 = h2_cc2aa(icc)
        idet = map_cc2aa(icc,3)
        !1 (p+p+|h+h+)
        tcc2(p1,p2,h1,h2,1) = +cic(idet,1)
        tcc2(p2,p1,h1,h2,1) = -cic(idet,1)
        tcc2(p1,p2,h2,h1,1) = -cic(idet,1)
        tcc2(p2,p1,h2,h1,1) = +cic(idet,1)
     end do
  end if
  if (nelact(2)>=2) then
     tcc2(:,:,:,:,2) = tcc2(:,:,:,:,1) !2 (p-p-|h-h-)
  end if

  if (nelact(2)>=1) then
     do icc = 1, ncc2ab
        p1 = p1_cc2ab(icc)
        p2 = p2_cc2ab(icc)
        h1 = h1_cc2ab(icc)
        h2 = h2_cc2ab(icc)
        idet = map_cc2ab(icc,3)     
        tcc2(p1,p2,h1,h2,3) = +cic(idet,1) !3 (p+p-|h+h-)     
        tcc2(p1,p2,h2,h1,5) = -cic(idet,1) !5 (p+p-|h-h+)
     end do
     tcc2(:,:,:,:,4) = tcc2(:,:,:,:,3) !4 (p-p+|h-h+)
     tcc2(:,:,:,:,6) = tcc2(:,:,:,:,5) !6 (p-p+|h+h-)
  end if
  !debug
  !do icc = 1, ncc2ab
  !   p1 = p1_cc2ab(icc)
  !   p2 = p2_cc2ab(icc)
  !   h1 = h1_cc2ab(icc)
  !   h2 = h2_cc2ab(icc)
  !   write(6,"('tdcc_gettcc2: ',8f20.10)") &
  !        dble(tcc2(p1,p2,h1,h2,3)), &
  !        dble(tcc2(p1,p2,h1,h2,4)), &
  !        dble(tcc2(p2,p1,h1,h2,5)), &
  !        dble(tcc2(p1,p2,h2,h1,6)), &
  !
  !        dble(tcc2(p2,p1,h2,h1,3)), &
  !        dble(tcc2(p2,p1,h2,h1,4)), &
  !        dble(tcc2(p1,p2,h2,h1,5)), &
  !        dble(tcc2(p2,p1,h1,h2,6))
  !end do
  !debug

end subroutine tdcc_gettcc2
!######################################################################
subroutine tdcc_getgcc2(cic, gcc2)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nact,ndetx,nelact
  use mod_cc, only : norb1,ncc2aa,ncc2ab,map_cc2aa,map_cc2ab
  use mod_cc, only : h1_cc2aa,h2_cc2aa,p1_cc2aa,p2_cc2aa
  use mod_cc, only : h1_cc2ab,h2_cc2ab,p1_cc2ab,p2_cc2ab
  implicit none
  complex(c_double_complex), intent(in) :: cic(1:ndetx,1:2)
  complex(c_double_complex), intent(out) :: &
       gcc2(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,1:6)
  integer(c_int) :: icc,p1,p2,h1,h2,idet

  gcc2 = 0d0
  if (nelact(3)< 2) return
  if (nelact(1)>=2) then
     do icc = 1, ncc2aa
        p1 = p1_cc2aa(icc)
        p2 = p2_cc2aa(icc)
        h1 = h1_cc2aa(icc)
        h2 = h2_cc2aa(icc)
        idet = map_cc2aa(icc,3)
        !1 (h+h+|p+p+)
        gcc2(h1,h2,p1,p2,1) = +cic(idet,2)
        gcc2(h2,h1,p1,p2,1) = -cic(idet,2)
        gcc2(h1,h2,p2,p1,1) = -cic(idet,2)
        gcc2(h2,h1,p2,p1,1) = +cic(idet,2)
     end do
  end if
  if (nelact(2)>=2) then
     gcc2(:,:,:,:,2) = gcc2(:,:,:,:,1) !2 (h-h-|p-p-)
  end if

  if (nelact(2)>=1) then
     do icc = 1, ncc2ab
        p1 = p1_cc2ab(icc)
        p2 = p2_cc2ab(icc)
        h1 = h1_cc2ab(icc)
        h2 = h2_cc2ab(icc)
        idet = map_cc2ab(icc,3)     
        gcc2(h1,h2,p1,p2,3) = +cic(idet,2) !3 (h+h-|p+p-)     
        gcc2(h1,h2,p2,p1,5) = -cic(idet,2) !5 (h+h-|p-p+)
     end do
     gcc2(:,:,:,:,4) = gcc2(:,:,:,:,3) !4 (h-h+|p-p+)
     gcc2(:,:,:,:,6) = gcc2(:,:,:,:,5) !6 (h-h+|p+p-)
  end if
  !debug
  !do icc = 1, ncc2ab
  !   p1 = p1_cc2ab(icc)
  !   p2 = p2_cc2ab(icc)
  !   h1 = h1_cc2ab(icc)
  !   h2 = h2_cc2ab(icc)
  !   write(6,"('tdcc_getgcc2: ',8f20.10)") &
  !        dble(gcc2(h1,h2,p1,p2,3)), &
  !        dble(gcc2(h1,h2,p1,p2,4)), &
  !        dble(gcc2(h1,h2,p2,p1,5)), &
  !        dble(gcc2(h2,h1,p1,p2,6)), &
! !
  !        dble(gcc2(h2,h1,p2,p1,3)), &
  !        dble(gcc2(h2,h1,p2,p1,4)), &
  !        dble(gcc2(h2,h1,p1,p2,5)), &
  !        dble(gcc2(h1,h2,p2,p1,6))
  !end do
  !debug

end subroutine tdcc_getgcc2
!######################################################################
subroutine tdcc_gettcc3(cic, tcc3)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nact,act1_ll,act1_ul,ndetx,nelact
  use mod_cc, only : norb1,ncc3aaa,ncc3aab,map_cc3aaa,map_cc3aab
  use mod_cc, only : h1_cc3aaa,h2_cc3aaa,h3_cc3aaa,p1_cc3aaa,p2_cc3aaa,p3_cc3aaa
  use mod_cc, only : h1_cc3aab,h2_cc3aab,h3_cc3aab,p1_cc3aab,p2_cc3aab,p3_cc3aab
  implicit none
  complex(c_double_complex), intent(in) :: cic(1:ndetx,1:2)
  complex(c_double_complex), intent(out) :: &
       tcc3((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,&
       1:norb1,1:norb1,1:norb1,1:20)
!  complex(c_double_complex), intent(out) :: &
!       tcc3((norb1+1):act1_ul,(norb1+1):act1_ul,(norb1+1):act1_ul,&
!       act1_ll:norb1,act1_ll:norb1,act1_ll:norb1,1:20)
  integer(c_int) :: icc,p1,p2,p3,h1,h2,h3,idet

  tcc3 = 0d0
  if (nelact(3)< 3) return
  if (nelact(1)>=3) then
     ! (+++|+++)
     do icc = 1, ncc3aaa
        p1 = p1_cc3aaa(icc)
        p2 = p2_cc3aaa(icc)
        p3 = p3_cc3aaa(icc)
        h1 = h1_cc3aaa(icc)
        h2 = h2_cc3aaa(icc)
        h3 = h3_cc3aaa(icc)
        idet = map_cc3aaa(icc,3)
        tcc3(p1,p2,p3,h1,h2,h3,1) = +cic(idet,1)
        tcc3(p1,p3,p2,h1,h2,h3,1) = -cic(idet,1)
        tcc3(p2,p1,p3,h1,h2,h3,1) = -cic(idet,1)
        tcc3(p3,p1,p2,h1,h2,h3,1) = +cic(idet,1)
        tcc3(p3,p2,p1,h1,h2,h3,1) = -cic(idet,1)
        tcc3(p2,p3,p1,h1,h2,h3,1) = +cic(idet,1)
     
        tcc3(p1,p2,p3,h1,h3,h2,1) = -cic(idet,1)
        tcc3(p1,p3,p2,h1,h3,h2,1) = +cic(idet,1)
        tcc3(p2,p1,p3,h1,h3,h2,1) = +cic(idet,1)
        tcc3(p3,p1,p2,h1,h3,h2,1) = -cic(idet,1)
        tcc3(p3,p2,p1,h1,h3,h2,1) = +cic(idet,1)
        tcc3(p2,p3,p1,h1,h3,h2,1) = -cic(idet,1)
     
        tcc3(p1,p2,p3,h2,h1,h3,1) = -cic(idet,1)
        tcc3(p1,p3,p2,h2,h1,h3,1) = +cic(idet,1)
        tcc3(p2,p1,p3,h2,h1,h3,1) = +cic(idet,1)
        tcc3(p3,p1,p2,h2,h1,h3,1) = -cic(idet,1)
        tcc3(p3,p2,p1,h2,h1,h3,1) = +cic(idet,1)
        tcc3(p2,p3,p1,h2,h1,h3,1) = -cic(idet,1)
     
        tcc3(p1,p2,p3,h3,h1,h2,1) = +cic(idet,1)
        tcc3(p1,p3,p2,h3,h1,h2,1) = -cic(idet,1)
        tcc3(p2,p1,p3,h3,h1,h2,1) = -cic(idet,1)
        tcc3(p3,p1,p2,h3,h1,h2,1) = +cic(idet,1)
        tcc3(p3,p2,p1,h3,h1,h2,1) = -cic(idet,1)
        tcc3(p2,p3,p1,h3,h1,h2,1) = +cic(idet,1)
     
        tcc3(p1,p2,p3,h3,h2,h1,1) = -cic(idet,1)
        tcc3(p1,p3,p2,h3,h2,h1,1) = +cic(idet,1)
        tcc3(p2,p1,p3,h3,h2,h1,1) = +cic(idet,1)
        tcc3(p3,p1,p2,h3,h2,h1,1) = -cic(idet,1)
        tcc3(p3,p2,p1,h3,h2,h1,1) = +cic(idet,1)
        tcc3(p2,p3,p1,h3,h2,h1,1) = -cic(idet,1)
     
        tcc3(p1,p2,p3,h2,h3,h1,1) = +cic(idet,1)
        tcc3(p1,p3,p2,h2,h3,h1,1) = -cic(idet,1)
        tcc3(p2,p1,p3,h2,h3,h1,1) = -cic(idet,1)
        tcc3(p3,p1,p2,h2,h3,h1,1) = +cic(idet,1)
        tcc3(p3,p2,p1,h2,h3,h1,1) = -cic(idet,1)
        tcc3(p2,p3,p1,h2,h3,h1,1) = +cic(idet,1)
     end do
  end if
  if (nelact(2)>=3) then
     ! (---|---)
     tcc3(:,:,:,:,:,:,2) = tcc3(:,:,:,:,:,:,1)
  end if
  if (nelact(1)>=2.and.nelact(2)>=1) then
     do icc = 1, ncc3aab
        p1 = p1_cc3aab(icc)
        p2 = p2_cc3aab(icc)
        p3 = p3_cc3aab(icc)
        h1 = h1_cc3aab(icc)
        h2 = h2_cc3aab(icc)
        h3 = h3_cc3aab(icc)
        idet = map_cc3aab(icc,3)
        ! (++-|++-)
        tcc3(p1,p2,p3,h1,h2,h3, 3) = +cic(idet,1)
        tcc3(p2,p1,p3,h1,h2,h3, 3) = -cic(idet,1)
        tcc3(p1,p2,p3,h2,h1,h3, 3) = -cic(idet,1)
        tcc3(p2,p1,p3,h2,h1,h3, 3) = +cic(idet,1)
     
        ! (++-|+-+)
        tcc3(p1,p2,p3,h1,h3,h2, 4) = -cic(idet,1)
        tcc3(p2,p1,p3,h1,h3,h2, 4) = +cic(idet,1)
        tcc3(p1,p2,p3,h2,h3,h1, 4) = +cic(idet,1)
        tcc3(p2,p1,p3,h2,h3,h1, 4) = -cic(idet,1)
     
        ! (++-|-++)
        tcc3(p1,p2,p3,h3,h2,h1, 5) = -cic(idet,1)
        tcc3(p2,p1,p3,h3,h2,h1, 5) = +cic(idet,1)
        tcc3(p1,p2,p3,h3,h1,h2, 5) = +cic(idet,1)
        tcc3(p2,p1,p3,h3,h1,h2, 5) = -cic(idet,1)
     
        ! (+-+|++-)
        tcc3(p1,p3,p2,h1,h2,h3, 6) = -cic(idet,1)
        tcc3(p2,p3,p1,h1,h2,h3, 6) = +cic(idet,1)
        tcc3(p1,p3,p2,h2,h1,h3, 6) = +cic(idet,1)
        tcc3(p2,p3,p1,h2,h1,h3, 6) = -cic(idet,1)
     
        ! (+-+|+-+)
        tcc3(p1,p3,p2,h1,h3,h2, 7) = +cic(idet,1)
        tcc3(p2,p3,p1,h1,h3,h2, 7) = -cic(idet,1)
        tcc3(p1,p3,p2,h2,h3,h1, 7) = -cic(idet,1)
        tcc3(p2,p3,p1,h2,h3,h1, 7) = +cic(idet,1)
     
        ! (+-+|-++)
        tcc3(p1,p3,p2,h3,h2,h1, 8) = +cic(idet,1)
        tcc3(p2,p3,p1,h3,h2,h1, 8) = -cic(idet,1)
        tcc3(p1,p3,p2,h3,h1,h2, 8) = -cic(idet,1)
        tcc3(p2,p3,p1,h3,h1,h2, 8) = +cic(idet,1)
     
        ! (-++|++-)
        tcc3(p3,p2,p1,h1,h2,h3, 9) = -cic(idet,1)
        tcc3(p3,p1,p2,h1,h2,h3, 9) = +cic(idet,1)
        tcc3(p3,p2,p1,h2,h1,h3, 9) = +cic(idet,1)
        tcc3(p3,p1,p2,h2,h1,h3, 9) = -cic(idet,1)
     
        ! (-++|+-+)
        tcc3(p3,p2,p1,h1,h3,h2,10) = +cic(idet,1)
        tcc3(p3,p1,p2,h1,h3,h2,10) = -cic(idet,1)
        tcc3(p3,p2,p1,h2,h3,h1,10) = -cic(idet,1)
        tcc3(p3,p1,p2,h2,h3,h1,10) = +cic(idet,1)
     
        ! (-++|-++)
        tcc3(p3,p2,p1,h3,h2,h1,11) = +cic(idet,1)
        tcc3(p3,p1,p2,h3,h2,h1,11) = -cic(idet,1)
        tcc3(p3,p2,p1,h3,h1,h2,11) = -cic(idet,1)
        tcc3(p3,p1,p2,h3,h1,h2,11) = +cic(idet,1)
     end do
  end if
  if (nelact(2)>=2) then
     tcc3(:,:,:,:,:,:,12) = tcc3(:,:,:,:,:,:, 3) ! (--+|--+)
     tcc3(:,:,:,:,:,:,13) = tcc3(:,:,:,:,:,:, 4) ! (--+|-+-)
     tcc3(:,:,:,:,:,:,14) = tcc3(:,:,:,:,:,:, 5) ! (--+|+--)
     
     tcc3(:,:,:,:,:,:,15) = tcc3(:,:,:,:,:,:, 6) ! (-+-|--+)
     tcc3(:,:,:,:,:,:,16) = tcc3(:,:,:,:,:,:, 7) ! (-+-|-+-)
     tcc3(:,:,:,:,:,:,17) = tcc3(:,:,:,:,:,:, 8) ! (-+-|+--)
     
     tcc3(:,:,:,:,:,:,18) = tcc3(:,:,:,:,:,:, 9) ! (+--|--+)
     tcc3(:,:,:,:,:,:,19) = tcc3(:,:,:,:,:,:,10) ! (+--|-+-)
     tcc3(:,:,:,:,:,:,20) = tcc3(:,:,:,:,:,:,11) ! (+--|+--)
  end if

end subroutine tdcc_gettcc3
!######################################################################
subroutine tdcc_getgcc3(cic, gcc3)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nact,act1_ll,act1_ul,ndetx
  use mod_cc, only : norb1,ncc3aaa,ncc3aab,map_cc3aaa,map_cc3aab
  use mod_cc, only : h1_cc3aaa,h2_cc3aaa,h3_cc3aaa,p1_cc3aaa,p2_cc3aaa,p3_cc3aaa
  use mod_cc, only : h1_cc3aab,h2_cc3aab,h3_cc3aab,p1_cc3aab,p2_cc3aab,p3_cc3aab
  implicit none
  complex(c_double_complex), intent(in) :: cic(1:ndetx,1:2)
  complex(c_double_complex), intent(out) :: &
       gcc3(1:norb1,1:norb1,1:norb1,&
       (norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:2)
!  complex(c_double_complex), intent(out) :: &
!       gcc3(act1_ll:norb1,act1_ll:norb1,act1_ll:norb1,&
!       (norb1+1):act1_ul,(norb1+1):act1_ul,(norb1+1):act1_ul,1:20)
  integer(c_int) :: icc,p1,p2,p3,h1,h2,h3,idet

  gcc3 = 0d0
  ! (+++|+++)
  do icc = 1, ncc3aaa
     p1 = p1_cc3aaa(icc)
     p2 = p2_cc3aaa(icc)
     p3 = p3_cc3aaa(icc)
     h1 = h1_cc3aaa(icc)
     h2 = h2_cc3aaa(icc)
     h3 = h3_cc3aaa(icc)
     idet = map_cc3aaa(icc,3)
     gcc3(h1,h2,h3,p1,p2,p3,1) = +cic(idet,2)
     gcc3(h1,h3,h2,p1,p2,p3,1) = -cic(idet,2)
     gcc3(h2,h1,h3,p1,p2,p3,1) = -cic(idet,2)
     gcc3(h3,h1,h2,p1,p2,p3,1) = +cic(idet,2)
     gcc3(h3,h2,h1,p1,p2,p3,1) = -cic(idet,2)
     gcc3(h2,h3,h1,p1,p2,p3,1) = +cic(idet,2)

     gcc3(h1,h2,h3,p1,p3,p2,1) = -cic(idet,2)
     gcc3(h1,h3,h2,p1,p3,p2,1) = +cic(idet,2)
     gcc3(h2,h1,h3,p1,p3,p2,1) = +cic(idet,2)
     gcc3(h3,h1,h2,p1,p3,p2,1) = -cic(idet,2)
     gcc3(h3,h2,h1,p1,p3,p2,1) = +cic(idet,2)
     gcc3(h2,h3,h1,p1,p3,p2,1) = -cic(idet,2)

     gcc3(h1,h2,h3,p2,p1,p3,1) = -cic(idet,2)
     gcc3(h1,h3,h2,p2,p1,p3,1) = +cic(idet,2)
     gcc3(h2,h1,h3,p2,p1,p3,1) = +cic(idet,2)
     gcc3(h3,h1,h2,p2,p1,p3,1) = -cic(idet,2)
     gcc3(h3,h2,h1,p2,p1,p3,1) = +cic(idet,2)
     gcc3(h2,h3,h1,p2,p1,p3,1) = -cic(idet,2)

     gcc3(h1,h2,h3,p3,p1,p2,1) = +cic(idet,2)
     gcc3(h1,h3,h2,p3,p1,p2,1) = -cic(idet,2)
     gcc3(h2,h1,h3,p3,p1,p2,1) = -cic(idet,2)
     gcc3(h3,h1,h2,p3,p1,p2,1) = +cic(idet,2)
     gcc3(h3,h2,h1,p3,p1,p2,1) = -cic(idet,2)
     gcc3(h2,h3,h1,p3,p1,p2,1) = +cic(idet,2)

     gcc3(h1,h2,h3,p3,p2,p1,1) = -cic(idet,2)
     gcc3(h1,h3,h2,p3,p2,p1,1) = +cic(idet,2)
     gcc3(h2,h1,h3,p3,p2,p1,1) = +cic(idet,2)
     gcc3(h3,h1,h2,p3,p2,p1,1) = -cic(idet,2)
     gcc3(h3,h2,h1,p3,p2,p1,1) = +cic(idet,2)
     gcc3(h2,h3,h1,p3,p2,p1,1) = -cic(idet,2)

     gcc3(h1,h2,h3,p2,p3,p1,1) = +cic(idet,2)
     gcc3(h1,h3,h2,p2,p3,p1,1) = -cic(idet,2)
     gcc3(h2,h1,h3,p2,p3,p1,1) = -cic(idet,2)
     gcc3(h3,h1,h2,p2,p3,p1,1) = +cic(idet,2)
     gcc3(h3,h2,h1,p2,p3,p1,1) = -cic(idet,2)
     gcc3(h2,h3,h1,p2,p3,p1,1) = +cic(idet,2)
  end do

  ! (---|---)
  gcc3(:,:,:,:,:,:,2) = gcc3(:,:,:,:,:,:,1)

  do icc = 1, ncc3aab
     p1 = p1_cc3aab(icc)
     p2 = p2_cc3aab(icc)
     p3 = p3_cc3aab(icc)
     h1 = h1_cc3aab(icc)
     h2 = h2_cc3aab(icc)
     h3 = h3_cc3aab(icc)
     idet = map_cc3aab(icc,3)
     ! (++-|++-)
     gcc3(h1,h2,h3,p1,p2,p3, 3) = +cic(idet,2)
     gcc3(h2,h1,h3,p1,p2,p3, 3) = -cic(idet,2)
     gcc3(h1,h2,h3,p2,p1,p3, 3) = -cic(idet,2)
     gcc3(h2,h1,h3,p2,p1,p3, 3) = +cic(idet,2)

     ! (++-|+-+)
     gcc3(h1,h2,h3,p1,p3,p2, 6) = -cic(idet,2)
     gcc3(h2,h1,h3,p1,p3,p2, 6) = +cic(idet,2)
     gcc3(h1,h2,h3,p2,p3,p1, 6) = +cic(idet,2)
     gcc3(h2,h1,h3,p2,p3,p1, 6) = -cic(idet,2)

     ! (++-|-++)
     gcc3(h1,h2,h3,p3,p2,p1, 9) = -cic(idet,2)
     gcc3(h2,h1,h3,p3,p2,p1, 9) = +cic(idet,2)
     gcc3(h1,h2,h3,p3,p1,p2, 9) = +cic(idet,2)
     gcc3(h2,h1,h3,p3,p1,p2, 9) = -cic(idet,2)

     ! (+-+|++-)
     gcc3(h1,h3,h2,p1,p2,p3, 4) = -cic(idet,2)
     gcc3(h2,h3,h1,p1,p2,p3, 4) = +cic(idet,2)
     gcc3(h1,h3,h2,p2,p1,p3, 4) = +cic(idet,2)
     gcc3(h2,h3,h1,p2,p1,p3, 4) = -cic(idet,2)

     ! (+-+|+-+)
     gcc3(h1,h3,h2,p1,p3,p2, 7) = +cic(idet,2)
     gcc3(h2,h3,h1,p1,p3,p2, 7) = -cic(idet,2)
     gcc3(h1,h3,h2,p2,p3,p1, 7) = -cic(idet,2)
     gcc3(h2,h3,h1,p2,p3,p1, 7) = +cic(idet,2)

     ! (+-+|-++)
     gcc3(h1,h3,h2,p3,p2,p1,10) = +cic(idet,2)
     gcc3(h2,h3,h1,p3,p2,p1,10) = -cic(idet,2)
     gcc3(h1,h3,h2,p3,p1,p2,10) = -cic(idet,2)
     gcc3(h2,h3,h1,p3,p1,p2,10) = +cic(idet,2)

     ! (-++|++-)
     gcc3(h3,h2,h1,p1,p2,p3, 5) = -cic(idet,2)
     gcc3(h3,h1,h2,p1,p2,p3, 5) = +cic(idet,2)
     gcc3(h3,h2,h1,p2,p1,p3, 5) = +cic(idet,2)
     gcc3(h3,h1,h2,p2,p1,p3, 5) = -cic(idet,2)

     ! (-++|+-+)
     gcc3(h3,h2,h1,p1,p3,p2, 8) = +cic(idet,2)
     gcc3(h3,h1,h2,p1,p3,p2, 8) = -cic(idet,2)
     gcc3(h3,h2,h1,p2,p3,p1, 8) = -cic(idet,2)
     gcc3(h3,h1,h2,p2,p3,p1, 8) = +cic(idet,2)

     ! (-++|-++)
     gcc3(h3,h2,h1,p3,p2,p1,11) = +cic(idet,2)
     gcc3(h3,h1,h2,p3,p2,p1,11) = -cic(idet,2)
     gcc3(h3,h2,h1,p3,p1,p2,11) = -cic(idet,2)
     gcc3(h3,h1,h2,p3,p1,p2,11) = +cic(idet,2)
  end do
  gcc3(:,:,:,:,:,:,12) = gcc3(:,:,:,:,:,:, 3) ! (--+|--+)
  gcc3(:,:,:,:,:,:,13) = gcc3(:,:,:,:,:,:, 4) ! (--+|-+-)
  gcc3(:,:,:,:,:,:,14) = gcc3(:,:,:,:,:,:, 5) ! (--+|+--)

  gcc3(:,:,:,:,:,:,15) = gcc3(:,:,:,:,:,:, 6) ! (-+-|--+)
  gcc3(:,:,:,:,:,:,16) = gcc3(:,:,:,:,:,:, 7) ! (-+-|-+-)
  gcc3(:,:,:,:,:,:,17) = gcc3(:,:,:,:,:,:, 8) ! (-+-|+--)

  gcc3(:,:,:,:,:,:,18) = gcc3(:,:,:,:,:,:, 9) ! (+--|--+)
  gcc3(:,:,:,:,:,:,19) = gcc3(:,:,:,:,:,:,10) ! (+--|-+-)
  gcc3(:,:,:,:,:,:,20) = gcc3(:,:,:,:,:,:,11) ! (+--|+--)

end subroutine tdcc_getgcc3
!######################################################################
subroutine tdcc_getgcc3_old(cic, gcc3)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nact,act1_ll,act1_ul,ndetx,nelact
  use mod_cc, only : norb1,ncc3aaa,ncc3aab,map_cc3aaa,map_cc3aab
  use mod_cc, only : h1_cc3aaa,h2_cc3aaa,h3_cc3aaa,p1_cc3aaa,p2_cc3aaa,p3_cc3aaa
  use mod_cc, only : h1_cc3aab,h2_cc3aab,h3_cc3aab,p1_cc3aab,p2_cc3aab,p3_cc3aab
  implicit none
  complex(c_double_complex), intent(in) :: cic(1:ndetx,1:2)
  complex(c_double_complex), intent(out) :: &
       gcc3(act1_ll:norb1,act1_ll:norb1,act1_ll:norb1,&
       (norb1+1):act1_ul,(norb1+1):act1_ul,(norb1+1):act1_ul,1:20)
  integer(c_int) :: icc,p1,p2,p3,h1,h2,h3,idet

  gcc3 = 0d0
  if (nelact(3)< 3) return
  if (nelact(1)>=3) then
     ! (+++|+++)
     do icc = 1, ncc3aaa
        p1 = p1_cc3aaa(icc)
        p2 = p2_cc3aaa(icc)
        p3 = p3_cc3aaa(icc)
        h1 = h1_cc3aaa(icc)
        h2 = h2_cc3aaa(icc)
        h3 = h3_cc3aaa(icc)
        idet = map_cc3aaa(icc,3)
        gcc3(h1,h2,h3,p1,p2,p3,1) = +cic(idet,2)
        gcc3(h1,h3,h2,p1,p2,p3,1) = -cic(idet,2)
        gcc3(h2,h1,h3,p1,p2,p3,1) = -cic(idet,2)
        gcc3(h3,h1,h2,p1,p2,p3,1) = +cic(idet,2)
        gcc3(h3,h2,h1,p1,p2,p3,1) = -cic(idet,2)
        gcc3(h2,h3,h1,p1,p2,p3,1) = +cic(idet,2)
     
        gcc3(h1,h2,h3,p1,p3,p2,1) = -cic(idet,2)
        gcc3(h1,h3,h2,p1,p3,p2,1) = +cic(idet,2)
        gcc3(h2,h1,h3,p1,p3,p2,1) = +cic(idet,2)
        gcc3(h3,h1,h2,p1,p3,p2,1) = -cic(idet,2)
        gcc3(h3,h2,h1,p1,p3,p2,1) = +cic(idet,2)
        gcc3(h2,h3,h1,p1,p3,p2,1) = -cic(idet,2)
     
        gcc3(h1,h2,h3,p2,p1,p3,1) = -cic(idet,2)
        gcc3(h1,h3,h2,p2,p1,p3,1) = +cic(idet,2)
        gcc3(h2,h1,h3,p2,p1,p3,1) = +cic(idet,2)
        gcc3(h3,h1,h2,p2,p1,p3,1) = -cic(idet,2)
        gcc3(h3,h2,h1,p2,p1,p3,1) = +cic(idet,2)
        gcc3(h2,h3,h1,p2,p1,p3,1) = -cic(idet,2)
     
        gcc3(h1,h2,h3,p3,p1,p2,1) = +cic(idet,2)
        gcc3(h1,h3,h2,p3,p1,p2,1) = -cic(idet,2)
        gcc3(h2,h1,h3,p3,p1,p2,1) = -cic(idet,2)
        gcc3(h3,h1,h2,p3,p1,p2,1) = +cic(idet,2)
        gcc3(h3,h2,h1,p3,p1,p2,1) = -cic(idet,2)
        gcc3(h2,h3,h1,p3,p1,p2,1) = +cic(idet,2)
     
        gcc3(h1,h2,h3,p3,p2,p1,1) = -cic(idet,2)
        gcc3(h1,h3,h2,p3,p2,p1,1) = +cic(idet,2)
        gcc3(h2,h1,h3,p3,p2,p1,1) = +cic(idet,2)
        gcc3(h3,h1,h2,p3,p2,p1,1) = -cic(idet,2)
        gcc3(h3,h2,h1,p3,p2,p1,1) = +cic(idet,2)
        gcc3(h2,h3,h1,p3,p2,p1,1) = -cic(idet,2)
     
        gcc3(h1,h2,h3,p2,p3,p1,1) = +cic(idet,2)
        gcc3(h1,h3,h2,p2,p3,p1,1) = -cic(idet,2)
        gcc3(h2,h1,h3,p2,p3,p1,1) = -cic(idet,2)
        gcc3(h3,h1,h2,p2,p3,p1,1) = +cic(idet,2)
        gcc3(h3,h2,h1,p2,p3,p1,1) = -cic(idet,2)
        gcc3(h2,h3,h1,p2,p3,p1,1) = +cic(idet,2)
     end do
  end if
  if (nelact(2)>=3) then
     ! (---|---)
     gcc3(:,:,:,:,:,:,2) = gcc3(:,:,:,:,:,:,1)
  end if
  if (nelact(1)>=2.and.nelact(2)>=1) then
     do icc = 1, ncc3aab
        p1 = p1_cc3aab(icc)
        p2 = p2_cc3aab(icc)
        p3 = p3_cc3aab(icc)
        h1 = h1_cc3aab(icc)
        h2 = h2_cc3aab(icc)
        h3 = h3_cc3aab(icc)
        idet = map_cc3aab(icc,3)
        ! (++-|++-)
        gcc3(h1,h2,h3,p1,p2,p3, 3) = +cic(idet,2)
        gcc3(h2,h1,h3,p1,p2,p3, 3) = -cic(idet,2)
        gcc3(h1,h2,h3,p2,p1,p3, 3) = -cic(idet,2)
        gcc3(h2,h1,h3,p2,p1,p3, 3) = +cic(idet,2)
     
        ! (++-|+-+)
        gcc3(h1,h2,h3,p1,p3,p2, 4) = -cic(idet,2)
        gcc3(h2,h1,h3,p1,p3,p2, 4) = +cic(idet,2)
        gcc3(h1,h2,h3,p2,p3,p1, 4) = +cic(idet,2)
        gcc3(h2,h1,h3,p2,p3,p1, 4) = -cic(idet,2)
     
        ! (++-|-++)
        gcc3(h1,h2,h3,p3,p2,p1, 5) = -cic(idet,2)
        gcc3(h2,h1,h3,p3,p2,p1, 5) = +cic(idet,2)
        gcc3(h1,h2,h3,p3,p1,p2, 5) = +cic(idet,2)
        gcc3(h2,h1,h3,p3,p1,p2, 5) = -cic(idet,2)
     
        ! (+-+|++-)
        gcc3(h1,h3,h2,p1,p2,p3, 6) = -cic(idet,2)
        gcc3(h2,h3,h1,p1,p2,p3, 6) = +cic(idet,2)
        gcc3(h1,h3,h2,p2,p1,p3, 6) = +cic(idet,2)
        gcc3(h2,h3,h1,p2,p1,p3, 6) = -cic(idet,2)
     
        ! (+-+|+-+)
        gcc3(h1,h3,h2,p1,p3,p2, 7) = +cic(idet,2)
        gcc3(h2,h3,h1,p1,p3,p2, 7) = -cic(idet,2)
        gcc3(h1,h3,h2,p2,p3,p1, 7) = -cic(idet,2)
        gcc3(h2,h3,h1,p2,p3,p1, 7) = +cic(idet,2)
     
        ! (+-+|-++)
        gcc3(h1,h3,h2,p3,p2,p1, 8) = +cic(idet,2)
        gcc3(h2,h3,h1,p3,p2,p1, 8) = -cic(idet,2)
        gcc3(h1,h3,h2,p3,p1,p2, 8) = -cic(idet,2)
        gcc3(h2,h3,h1,p3,p1,p2, 8) = +cic(idet,2)
     
        ! (-++|++-)
        gcc3(h3,h2,h1,p1,p2,p3, 9) = -cic(idet,2)
        gcc3(h3,h1,h2,p1,p2,p3, 9) = +cic(idet,2)
        gcc3(h3,h2,h1,p2,p1,p3, 9) = +cic(idet,2)
        gcc3(h3,h1,h2,p2,p1,p3, 9) = -cic(idet,2)
     
        ! (-++|+-+)
        gcc3(h3,h2,h1,p1,p3,p2,10) = +cic(idet,2)
        gcc3(h3,h1,h2,p1,p3,p2,10) = -cic(idet,2)
        gcc3(h3,h2,h1,p2,p3,p1,10) = -cic(idet,2)
        gcc3(h3,h1,h2,p2,p3,p1,10) = +cic(idet,2)
     
        ! (-++|-++)
        gcc3(h3,h2,h1,p3,p2,p1,11) = +cic(idet,2)
        gcc3(h3,h1,h2,p3,p2,p1,11) = -cic(idet,2)
        gcc3(h3,h2,h1,p3,p1,p2,11) = -cic(idet,2)
        gcc3(h3,h1,h2,p3,p1,p2,11) = +cic(idet,2)
     end do
  end if
  if (nelact(2)>=2) then
     gcc3(:,:,:,:,:,:,12) = gcc3(:,:,:,:,:,:, 3) ! (--+|--+)
     gcc3(:,:,:,:,:,:,13) = gcc3(:,:,:,:,:,:, 4) ! (--+|-+-)
     gcc3(:,:,:,:,:,:,14) = gcc3(:,:,:,:,:,:, 5) ! (--+|+--)
     
     gcc3(:,:,:,:,:,:,15) = gcc3(:,:,:,:,:,:, 6) ! (-+-|--+)
     gcc3(:,:,:,:,:,:,16) = gcc3(:,:,:,:,:,:, 7) ! (-+-|-+-)
     gcc3(:,:,:,:,:,:,17) = gcc3(:,:,:,:,:,:, 8) ! (-+-|+--)
     
     gcc3(:,:,:,:,:,:,18) = gcc3(:,:,:,:,:,:, 9) ! (+--|--+)
     gcc3(:,:,:,:,:,:,19) = gcc3(:,:,:,:,:,:,10) ! (+--|-+-)
     gcc3(:,:,:,:,:,:,20) = gcc3(:,:,:,:,:,:,11) ! (+--|+--)
  end if

end subroutine tdcc_getgcc3_old
!######################################################################
