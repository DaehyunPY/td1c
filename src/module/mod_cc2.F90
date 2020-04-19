!################################################################################
module mod_cc2

  use, intrinsic :: iso_c_binding
  implicit none

  public

  integer(c_int) :: spin_focka
  integer(c_int) :: spin_int2aa
  integer(c_int) :: spin_int2ab
  integer(c_int) :: spin_t2aa
  integer(c_int) :: spin_t2ab
  integer(c_int) :: spin_t3aaa
  integer(c_int) :: spin_t3aab
  integer(c_int) :: spin_g1a
  integer(c_int) :: spin_g2aa
  integer(c_int) :: spin_g2ab
  integer(c_int) :: spin_g3aaa
  integer(c_int) :: spin_g3aab
  complex(c_double_complex), allocatable :: cc_work1(:)
  complex(c_double_complex), allocatable :: cc_work2(:)
  complex(c_double_complex), allocatable :: cc_work3(:)

  integer(c_int) :: noo,nov,nvv
  integer(c_int),allocatable :: h1_oo(:),h2_oo(:)
  integer(c_int),allocatable :: h1_ov(:),p2_ov(:)
  integer(c_int),allocatable :: p1_vv(:),p2_vv(:)
  integer(c_int) :: nooooaa,nooovaa,noovvaa,novovaa,novvvaa,nvvvvaa
  integer(c_int),allocatable :: h1_ooooaa(:),h2_ooooaa(:),h3_ooooaa(:),h4_ooooaa(:)
  integer(c_int),allocatable :: h1_ooovaa(:),h2_ooovaa(:),h3_ooovaa(:),p4_ooovaa(:)
  integer(c_int),allocatable :: h1_oovvaa(:),h2_oovvaa(:),p3_oovvaa(:),p4_oovvaa(:)
  integer(c_int),allocatable :: h1_ovovaa(:),p2_ovovaa(:),h3_ovovaa(:),p4_ovovaa(:)
  integer(c_int),allocatable :: h1_ovvvaa(:),p2_ovvvaa(:),p3_ovvvaa(:),p4_ovvvaa(:)
  integer(c_int),allocatable :: p1_vvvvaa(:),p2_vvvvaa(:),p3_vvvvaa(:),p4_vvvvaa(:)
  integer(c_int) :: nooooab,nooovab,noovvab,novovab,novvvab,nvvvvab
  integer(c_int),allocatable :: h1_ooooab(:),h2_ooooab(:),h3_ooooab(:),h4_ooooab(:)
  integer(c_int),allocatable :: h1_ooovab(:),h2_ooovab(:),h3_ooovab(:),p4_ooovab(:)
  integer(c_int),allocatable :: h1_oovvab(:),h2_oovvab(:),p3_oovvab(:),p4_oovvab(:)
  integer(c_int),allocatable :: h1_ovovab(:),p2_ovovab(:),h3_ovovab(:),p4_ovovab(:)
  integer(c_int),allocatable :: h1_ovvvab(:),p2_ovvvab(:),p3_ovvvab(:),p4_ovvvab(:)
  integer(c_int),allocatable :: p1_vvvvab(:),p2_vvvvab(:),p3_vvvvab(:),p4_vvvvab(:)

  contains
  subroutine mod_cc2_index_break
!nyi    icc = -1
!nyi    a = -1
!nyi    b = -1
!nyi    c = -1
!nyi    d = -1
!nyi    e = -1
!nyi    i = -1
!nyi    j = -1
!nyi    k = -1
!nyi    l = -1
!nyi    m = -1
  end subroutine mod_cc2_index_break
end module mod_cc2
!################################################################################
