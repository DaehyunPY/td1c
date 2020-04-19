!######################################################################
subroutine tdcc_print(cic)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nact
  use mod_cc, only : cc_rank
  implicit none
  complex(c_double_complex), intent(in) :: cic(1:*)

  call tdcc_printcc0(cic)
  if (cc_rank >= 1) call tdcc_printcc1(cic)
  if (cc_rank >= 2) call tdcc_printcc2(cic)
  if (cc_rank >= 3) call tdcc_printcc3(cic)

end subroutine tdcc_print
!######################################################################
subroutine tdcc_printcc0(cic)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : ndetx
  use mod_cc, only : map_cc0
  implicit none
  complex(c_double_complex), intent(in) :: cic(1:ndetx,1:2)

  write(6,"('tdcc_0: ',4f20.10)") &
       cic(map_cc0(1,3),1),cic(map_cc0(1,3),2)

end subroutine tdcc_printcc0
!######################################################################
subroutine tdcc_printcc1(cic)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nact,ndetx
  use mod_cc, only : norb1,ncc1a
  use mod_cc, only : h1_cc1a,p1_cc1a,map_cc1a
  implicit none
  complex(c_double_complex), intent(in) :: cic(1:ndetx,1:2)
  integer(c_long) :: icc,p1,h1,idet

  do icc = 1, ncc1a
     p1 = p1_cc1a(icc)
     h1 = h1_cc1a(icc)
     idet = map_cc1a(icc,3)
     write(6,"('tdcc_1aa: ',2i5,4f20.10)") p1,h1,&
          cic(idet,1),cic(idet,2)
  end do

end subroutine tdcc_printcc1
!######################################################################
subroutine tdcc_printcc2(cic)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nact,ndetx
  use mod_cc, only : norb1,ncc2aa,ncc2ab,map_cc2aa,map_cc2ab
  use mod_cc, only : h1_cc2aa,h2_cc2aa,p1_cc2aa,p2_cc2aa
  use mod_cc, only : h1_cc2ab,h2_cc2ab,p1_cc2ab,p2_cc2ab
  implicit none
  complex(c_double_complex), intent(in) :: cic(1:ndetx,1:2)
  integer(c_long) :: icc,p1,p2,h1,h2,idet

  do icc = 1, ncc2aa
     p1 = p1_cc2aa(icc)
     p2 = p2_cc2aa(icc)
     h1 = h1_cc2aa(icc)
     h2 = h2_cc2aa(icc)
     idet = map_cc2aa(icc,3)
     write(6,"('tdcc_2aaaa: ',4i5,4f20.10)") p1,p2,h1,h2,&
          cic(idet,1),cic(idet,2)
  end do
  do icc = 1, ncc2ab
     p1 = p1_cc2ab(icc)
     p2 = p2_cc2ab(icc)
     h1 = h1_cc2ab(icc)
     h2 = h2_cc2ab(icc)
     idet = map_cc2ab(icc,3)
     write(6,"('tdcc_2abab: ',4i5,4f20.10)") p1,p2,h1,h2,&
          cic(idet,1),cic(idet,2)
  end do

end subroutine tdcc_printcc2
!######################################################################
subroutine tdcc_printcc3(cic)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nact,ndetx
  use mod_cc, only : norb1,ncc3aaa,ncc3aab,map_cc3aaa,map_cc3aab
  use mod_cc, only : h1_cc3aaa,h2_cc3aaa,h3_cc3aaa,p1_cc3aaa,p2_cc3aaa,p3_cc3aaa
  use mod_cc, only : h1_cc3aab,h2_cc3aab,h3_cc3aab,p1_cc3aab,p2_cc3aab,p3_cc3aab
  implicit none
  complex(c_double_complex), intent(in) :: cic(1:ndetx,1:2)
  integer(c_long) :: icc,p1,p2,p3,h1,h2,h3,idet

  do icc = 1, ncc3aaa
     p1 = p1_cc3aaa(icc)
     p2 = p2_cc3aaa(icc)
     p3 = p3_cc3aaa(icc)
     h1 = h1_cc3aaa(icc)
     h2 = h2_cc3aaa(icc)
     h3 = h3_cc3aaa(icc)
     idet = map_cc3aaa(icc,3)
     write(6,"('tdcc_3aaaaaa: ',6i5,4f20.10)") p1,p2,p3,h1,h2,h3,&
          cic(idet,1),cic(idet,2)
  end do

  do icc = 1, ncc3aab
     p1 = p1_cc3aab(icc)
     p2 = p2_cc3aab(icc)
     p3 = p3_cc3aab(icc)
     h1 = h1_cc3aab(icc)
     h2 = h2_cc3aab(icc)
     h3 = h3_cc3aab(icc)
     idet = map_cc3aab(icc,3)
     write(6,"('tdcc_3aabaab: ',6i5,4f20.10)") p1,p2,p3,h1,h2,h3,&
          cic(idet,1),cic(idet,2)
  end do

end subroutine tdcc_printcc3
!######################################################################
