!######################################################################
subroutine tdcc_putcc(cic)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nact
  use mod_cc, only : cc_rank,t0out,g0out,t1out,g1out,t2out,g2out,t3out,g3out
  implicit none
  complex(c_double_complex), intent(out) :: cic(1:*)

  call tdcc_puttcc0(t0out,cic)
  if (cc_rank >= 1) call tdcc_puttcc1(t1out,cic)
  if (cc_rank >= 2) call tdcc_puttcc2(t2out,cic)
  if (cc_rank >= 3) call tdcc_puttcc3(t3out,cic)

  call tdcc_putgcc0(g0out,cic)
  if (cc_rank >= 1) call tdcc_putgcc1(g1out,cic)
  if (cc_rank >= 2) call tdcc_putgcc2(g2out,cic)
  if (cc_rank >= 3) call tdcc_putgcc3(g3out,cic)

end subroutine tdcc_putcc
!######################################################################
subroutine tdcc_puttcc0(tcc0, cic)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : ndetx
  use mod_cc, only : map_cc0
  implicit none
  complex(c_double_complex), intent(in) :: tcc0
  complex(c_double_complex), intent(out) :: cic(1:ndetx,1:2)

  cic(map_cc0(1,3),1) = tcc0

end subroutine tdcc_puttcc0
!######################################################################
subroutine tdcc_putgcc0(gcc0, cic)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : ndetx
  use mod_cc, only : map_cc0
  implicit none
  complex(c_double_complex), intent(in) :: gcc0
  complex(c_double_complex), intent(out) :: cic(1:ndetx,1:2)

  cic(map_cc0(1,3),2) = gcc0

end subroutine tdcc_putgcc0
!######################################################################
subroutine tdcc_puttcc1(tcc1, cic)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nact,ndetx
  use mod_cc, only : norb1,ncc1a
  use mod_cc, only : h1_cc1a,p1_cc1a,map_cc1a
  implicit none
  complex(c_double_complex), intent(in) :: tcc1((norb1+1):nact,1:norb1,1:2)
  complex(c_double_complex), intent(out) :: cic(1:ndetx,1:2)
  integer(c_int) :: icc,p1,h1,idet

  do icc = 1, ncc1a
     p1 = p1_cc1a(icc)
     h1 = h1_cc1a(icc)
     idet = map_cc1a(icc,3)
     cic(idet,1) = tcc1(p1,h1,1)
  end do

end subroutine tdcc_puttcc1
!######################################################################
subroutine tdcc_putgcc1(gcc1, cic)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nact,ndetx
  use mod_cc, only : norb1,ncc1a
  use mod_cc, only : h1_cc1a,p1_cc1a,map_cc1a
  implicit none
  complex(c_double_complex), intent(in) :: gcc1(1:norb1,(norb1+1):nact,1:2)
  complex(c_double_complex), intent(out) :: cic(1:ndetx,1:2)
  integer(c_int) :: icc,p1,h1,idet

  do icc = 1, ncc1a
     p1 = p1_cc1a(icc)
     h1 = h1_cc1a(icc)
     idet = map_cc1a(icc,3)
     cic(idet,2) = gcc1(h1,p1,1)
  end do

end subroutine tdcc_putgcc1
!######################################################################
subroutine tdcc_puttcc2(tcc2, cic)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nact,ndetx
  use mod_cc, only : norb1,ncc2aa,ncc2ab,map_cc2aa,map_cc2ab
  use mod_cc, only : h1_cc2aa,h2_cc2aa,p1_cc2aa,p2_cc2aa
  use mod_cc, only : h1_cc2ab,h2_cc2ab,p1_cc2ab,p2_cc2ab
  implicit none
  complex(c_double_complex), intent(in) :: tcc2((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:2)
  complex(c_double_complex), intent(out) :: cic(1:ndetx,1:2)
  integer(c_int) :: icc,p1,p2,h1,h2,idet

  ! (p+p+|h+h+)
  do icc = 1, ncc2aa
     p1 = p1_cc2aa(icc)
     p2 = p2_cc2aa(icc)
     h1 = h1_cc2aa(icc)
     h2 = h2_cc2aa(icc)
     idet = map_cc2aa(icc,3)
     cic(idet,1) = tcc2(p1,p2,h1,h2,1)
  end do

  ! (p+p-|h+h-)
  do icc = 1, ncc2ab
     p1 = p1_cc2ab(icc)
     p2 = p2_cc2ab(icc)
     h1 = h1_cc2ab(icc)
     h2 = h2_cc2ab(icc)
     idet = map_cc2ab(icc,3)
     cic(idet,1) = tcc2(p1,p2,h1,h2,2)
  end do

end subroutine tdcc_puttcc2
!######################################################################
subroutine tdcc_putgcc2(gcc2, cic)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nact,ndetx
  use mod_cc, only : norb1,ncc2aa,ncc2ab,map_cc2aa,map_cc2ab
  use mod_cc, only : h1_cc2aa,h2_cc2aa,p1_cc2aa,p2_cc2aa
  use mod_cc, only : h1_cc2ab,h2_cc2ab,p1_cc2ab,p2_cc2ab
  implicit none
  complex(c_double_complex), intent(in) :: gcc2(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,1:2)
  complex(c_double_complex), intent(out) :: cic(1:ndetx,1:2)
  integer(c_int) :: icc,p1,p2,h1,h2,idet

  ! (h+h+|p+p+)
  do icc = 1, ncc2aa
     p1 = p1_cc2aa(icc)
     p2 = p2_cc2aa(icc)
     h1 = h1_cc2aa(icc)
     h2 = h2_cc2aa(icc)
     idet = map_cc2aa(icc,3)
     cic(idet,2) = gcc2(h1,h2,p1,p2,1)
  end do

  ! (h+h-|p+p-)
  do icc = 1, ncc2ab
     p1 = p1_cc2ab(icc)
     p2 = p2_cc2ab(icc)
     h1 = h1_cc2ab(icc)
     h2 = h2_cc2ab(icc)
     idet = map_cc2ab(icc,3)
     cic(idet,2) = gcc2(h1,h2,p1,p2,2)
  end do

end subroutine tdcc_putgcc2
!######################################################################
subroutine tdcc_puttcc3(tcc3, cic)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nact,act1_ll,act1_ul,ndetx
  use mod_cc, only : norb1,ncc3aaa,ncc3aab,map_cc3aaa,map_cc3aab
  use mod_cc, only : h1_cc3aaa,h2_cc3aaa,h3_cc3aaa,p1_cc3aaa,p2_cc3aaa,p3_cc3aaa
  use mod_cc, only : h1_cc3aab,h2_cc3aab,h3_cc3aab,p1_cc3aab,p2_cc3aab,p3_cc3aab
  implicit none
  complex(c_double_complex), intent(in) :: &
       tcc3((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,&
       1:norb1,1:norb1,1:norb1,1:20)
!  complex(c_double_complex), intent(in) :: &
!       tcc3((norb1+1):act1_ul,(norb1+1):act1_ul,&
!       (norb1+1):act1_ul,act1_ll:norb1,act1_ll:norb1,act1_ll:norb1,1:2)
  complex(c_double_complex), intent(out) :: cic(1:ndetx,1:2)
  integer(c_int) :: icc,p1,p2,p3,h1,h2,h3,idet

  ! (+++|+++)
  do icc = 1, ncc3aaa
     p1 = p1_cc3aaa(icc)
     p2 = p2_cc3aaa(icc)
     p3 = p3_cc3aaa(icc)
     h1 = h1_cc3aaa(icc)
     h2 = h2_cc3aaa(icc)
     h3 = h3_cc3aaa(icc)
     idet = map_cc3aaa(icc,3)
     cic(idet,1) = tcc3(p1,p2,p3,h1,h2,h3,1)
  end do

  ! (++-|++-)
  do icc = 1, ncc3aab
     p1 = p1_cc3aab(icc)
     p2 = p2_cc3aab(icc)
     p3 = p3_cc3aab(icc)
     h1 = h1_cc3aab(icc)
     h2 = h2_cc3aab(icc)
     h3 = h3_cc3aab(icc)
     idet = map_cc3aab(icc,3)
     cic(idet,1) = tcc3(p1,p2,p3,h1,h2,h3,2)
  end do

end subroutine tdcc_puttcc3
!######################################################################
subroutine tdcc_putgcc3(gcc3, cic)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nact,act1_ll,act1_ul,ndetx
  use mod_cc, only : norb1,ncc3aaa,ncc3aab,map_cc3aaa,map_cc3aab
  use mod_cc, only : h1_cc3aaa,h2_cc3aaa,h3_cc3aaa,p1_cc3aaa,p2_cc3aaa,p3_cc3aaa
  use mod_cc, only : h1_cc3aab,h2_cc3aab,h3_cc3aab,p1_cc3aab,p2_cc3aab,p3_cc3aab
  implicit none
  complex(c_double_complex), intent(in) :: &
       gcc3(1:norb1,1:norb1,1:norb1,&
       (norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:2)
!  complex(c_double_complex), intent(in) :: &
!       gcc3(act1_ll:norb1,act1_ll:norb1,act1_ll:norb1,&
!       (norb1+1):act1_ul,(norb1+1):act1_ul,(norb1+1):act1_ul,1:2)
  complex(c_double_complex), intent(out) :: cic(1:ndetx,1:2)
  integer(c_int) :: icc,p1,p2,p3,h1,h2,h3,idet

  ! (+++|+++)
  do icc = 1, ncc3aaa
     p1 = p1_cc3aaa(icc)
     p2 = p2_cc3aaa(icc)
     p3 = p3_cc3aaa(icc)
     h1 = h1_cc3aaa(icc)
     h2 = h2_cc3aaa(icc)
     h3 = h3_cc3aaa(icc)
     idet = map_cc3aaa(icc,3)
     cic(idet,2) = gcc3(h1,h2,h3,p1,p2,p3,1)
  end do

  ! (++-|++-)
  do icc = 1, ncc3aab
     p1 = p1_cc3aab(icc)
     p2 = p2_cc3aab(icc)
     p3 = p3_cc3aab(icc)
     h1 = h1_cc3aab(icc)
     h2 = h2_cc3aab(icc)
     h3 = h3_cc3aab(icc)
     idet = map_cc3aab(icc,3)
     cic(idet,2) = gcc3(h1,h2,h3,p1,p2,p3,2)
  end do

end subroutine tdcc_putgcc3
!######################################################################
