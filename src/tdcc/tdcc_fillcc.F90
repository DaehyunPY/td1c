!######################################################################
subroutine tdcc_fillooooaa(i0)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nact
  use mod_cc, only : norb1
  use mod_cc2

  implicit none
  complex(c_double_complex), intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,1:norb1)
  integer(c_long) :: icc,i,j,k,l

  do icc = 1, nooooaa
     i = h1_ooooaa(icc)
     j = h2_ooooaa(icc)
     k = h3_ooooaa(icc)
     l = h4_ooooaa(icc)
     i0(j,i,k,l) = -i0(i,j,k,l)
     i0(i,j,l,k) = -i0(i,j,k,l)
     i0(j,i,l,k) =  i0(i,j,k,l)
  end do

end subroutine tdcc_fillooooaa
!######################################################################
subroutine tdcc_fillooovaa(i0)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nact
  use mod_cc, only : norb1
  use mod_cc2

  implicit none
  complex(c_double_complex), intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer(c_long) :: icc,i,j,k,a

  do icc = 1, nooovaa
     i = h1_ooovaa(icc)
     j = h2_ooovaa(icc)
     k = h3_ooovaa(icc)
     a = p4_ooovaa(icc)
     i0(j,i,k,a) = -i0(i,j,k,a)
  end do

end subroutine tdcc_fillooovaa
!######################################################################
subroutine tdcc_filloovoaa(i0)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nact
  use mod_cc, only : norb1
  use mod_cc2

  implicit none
  complex(c_double_complex), intent(inout) :: &
       i0(1:norb1,1:norb1,(norb1+1):nact,1:norb1)       
  integer(c_long) :: icc,i,j,a,k

  do icc = 1, nooovaa
     i = h1_ooovaa(icc)
     j = h2_ooovaa(icc)
     a = p4_ooovaa(icc)
     k = h3_ooovaa(icc)
     i0(j,i,a,k) = -i0(i,j,a,k)
  end do

end subroutine tdcc_filloovoaa
!######################################################################
subroutine tdcc_fillovooaa(i0)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nact
  use mod_cc, only : norb1
  use mod_cc2

  implicit none
  complex(c_double_complex), intent(inout) :: &
       i0(1:norb1,(norb1+1):nact,1:norb1,1:norb1)       
  integer(c_long) :: icc,i,j,a,k

  do icc = 1, nooovaa
     i = h3_ooovaa(icc)
     a = p4_ooovaa(icc)
     j = h1_ooovaa(icc)
     k = h2_ooovaa(icc)
     i0(i,a,k,j) = -i0(i,a,j,k)
  end do

end subroutine tdcc_fillovooaa
!######################################################################
subroutine tdcc_fillvvovaa(i0)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nact
  use mod_cc, only : norb1
  use mod_cc2

  implicit none
  complex(c_double_complex), intent(inout) :: &
       i0((norb1+1):nact,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer(c_long) :: icc,i,a,b,c

  do icc = 1, novvvaa
     a = p3_ovvvaa(icc)
     b = p4_ovvvaa(icc)
     i = h1_ovvvaa(icc)
     c = p2_ovvvaa(icc)
     i0(b,a,i,c) = -i0(a,b,i,c)
  end do

end subroutine tdcc_fillvvovaa
!######################################################################
subroutine tdcc_fillovvvaa(i0)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nact
  use mod_cc, only : norb1
  use mod_cc2

  implicit none
  complex(c_double_complex), intent(inout) :: &
       i0(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  integer(c_long) :: icc,i,a,b,c

  do icc = 1, novvvaa
     i = h1_ovvvaa(icc)
     a = p2_ovvvaa(icc)
     b = p3_ovvvaa(icc)
     c = p4_ovvvaa(icc)
     i0(i,a,c,b) = -i0(i,a,b,c)
  end do

end subroutine tdcc_fillovvvaa
!######################################################################
subroutine tdcc_fillvvvvaa(i0)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nact
  use mod_cc, only : norb1
  use mod_cc2

  implicit none
  complex(c_double_complex), intent(inout) :: &
       i0((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  integer(c_long) :: icc,a,b,c,d

  do icc = 1, nvvvvaa
     a = p1_vvvvaa(icc)
     b = p2_vvvvaa(icc)
     c = p3_vvvvaa(icc)
     d = p4_vvvvaa(icc)
     i0(b,a,c,d) = -i0(a,b,c,d)
     i0(a,b,d,c) = -i0(a,b,c,d)
     i0(b,a,d,c) =  i0(a,b,c,d)
  end do

end subroutine tdcc_fillvvvvaa
!######################################################################
subroutine tdcc_filltcc3aaa(tcc3)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nact
  use mod_cc, only : norb1,ncc3aaa
  use mod_cc, only : h1_cc3aaa,h2_cc3aaa,h3_cc3aaa,p1_cc3aaa,p2_cc3aaa,p3_cc3aaa
  implicit none
  complex(c_double_complex), intent(inout) :: &
       tcc3((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)
  integer(c_long) :: icc,p1,p2,p3,h1,h2,h3

  ! (+++|+++)
  do icc = 1, ncc3aaa
     p1 = p1_cc3aaa(icc)
     p2 = p2_cc3aaa(icc)
     p3 = p3_cc3aaa(icc)
     h1 = h1_cc3aaa(icc)
     h2 = h2_cc3aaa(icc)
     h3 = h3_cc3aaa(icc)
     !tcc3(p1,p2,p3,h1,h2,h3) = +cic(idet)
     tcc3(p1,p3,p2,h1,h2,h3) = -tcc3(p1,p2,p3,h1,h2,h3)
     tcc3(p2,p1,p3,h1,h2,h3) = -tcc3(p1,p2,p3,h1,h2,h3)
     tcc3(p3,p1,p2,h1,h2,h3) = +tcc3(p1,p2,p3,h1,h2,h3)
     tcc3(p3,p2,p1,h1,h2,h3) = -tcc3(p1,p2,p3,h1,h2,h3)
     tcc3(p2,p3,p1,h1,h2,h3) = +tcc3(p1,p2,p3,h1,h2,h3)

     tcc3(p1,p2,p3,h1,h3,h2) = -tcc3(p1,p2,p3,h1,h2,h3)
     tcc3(p1,p3,p2,h1,h3,h2) = +tcc3(p1,p2,p3,h1,h2,h3)
     tcc3(p2,p1,p3,h1,h3,h2) = +tcc3(p1,p2,p3,h1,h2,h3)
     tcc3(p3,p1,p2,h1,h3,h2) = -tcc3(p1,p2,p3,h1,h2,h3)
     tcc3(p3,p2,p1,h1,h3,h2) = +tcc3(p1,p2,p3,h1,h2,h3)
     tcc3(p2,p3,p1,h1,h3,h2) = -tcc3(p1,p2,p3,h1,h2,h3)

     tcc3(p1,p2,p3,h2,h1,h3) = -tcc3(p1,p2,p3,h1,h2,h3)
     tcc3(p1,p3,p2,h2,h1,h3) = +tcc3(p1,p2,p3,h1,h2,h3)
     tcc3(p2,p1,p3,h2,h1,h3) = +tcc3(p1,p2,p3,h1,h2,h3)
     tcc3(p3,p1,p2,h2,h1,h3) = -tcc3(p1,p2,p3,h1,h2,h3)
     tcc3(p3,p2,p1,h2,h1,h3) = +tcc3(p1,p2,p3,h1,h2,h3)
     tcc3(p2,p3,p1,h2,h1,h3) = -tcc3(p1,p2,p3,h1,h2,h3)

     tcc3(p1,p2,p3,h3,h1,h2) = +tcc3(p1,p2,p3,h1,h2,h3)
     tcc3(p1,p3,p2,h3,h1,h2) = -tcc3(p1,p2,p3,h1,h2,h3)
     tcc3(p2,p1,p3,h3,h1,h2) = -tcc3(p1,p2,p3,h1,h2,h3)
     tcc3(p3,p1,p2,h3,h1,h2) = +tcc3(p1,p2,p3,h1,h2,h3)
     tcc3(p3,p2,p1,h3,h1,h2) = -tcc3(p1,p2,p3,h1,h2,h3)
     tcc3(p2,p3,p1,h3,h1,h2) = +tcc3(p1,p2,p3,h1,h2,h3)

     tcc3(p1,p2,p3,h3,h2,h1) = -tcc3(p1,p2,p3,h1,h2,h3)
     tcc3(p1,p3,p2,h3,h2,h1) = +tcc3(p1,p2,p3,h1,h2,h3)
     tcc3(p2,p1,p3,h3,h2,h1) = +tcc3(p1,p2,p3,h1,h2,h3)
     tcc3(p3,p1,p2,h3,h2,h1) = -tcc3(p1,p2,p3,h1,h2,h3)
     tcc3(p3,p2,p1,h3,h2,h1) = +tcc3(p1,p2,p3,h1,h2,h3)
     tcc3(p2,p3,p1,h3,h2,h1) = -tcc3(p1,p2,p3,h1,h2,h3)

     tcc3(p1,p2,p3,h2,h3,h1) = +tcc3(p1,p2,p3,h1,h2,h3)
     tcc3(p1,p3,p2,h2,h3,h1) = -tcc3(p1,p2,p3,h1,h2,h3)
     tcc3(p2,p1,p3,h2,h3,h1) = -tcc3(p1,p2,p3,h1,h2,h3)
     tcc3(p3,p1,p2,h2,h3,h1) = +tcc3(p1,p2,p3,h1,h2,h3)
     tcc3(p3,p2,p1,h2,h3,h1) = -tcc3(p1,p2,p3,h1,h2,h3)
     tcc3(p2,p3,p1,h2,h3,h1) = +tcc3(p1,p2,p3,h1,h2,h3)
  end do

end subroutine tdcc_filltcc3aaa
!######################################################################
subroutine tdcc_filltcc3aab(tcc3)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nact
  use mod_cc, only : norb1,ncc3aab
  use mod_cc, only : h1_cc3aab,h2_cc3aab,h3_cc3aab,p1_cc3aab,p2_cc3aab,p3_cc3aab
  implicit none
  complex(c_double_complex), intent(inout) :: &
       tcc3((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)
  integer(c_long) :: icc,p1,p2,p3,h1,h2,h3,idet

  do icc = 1, ncc3aab
     p1 = p1_cc3aab(icc)
     p2 = p2_cc3aab(icc)
     p3 = p3_cc3aab(icc)
     h1 = h1_cc3aab(icc)
     h2 = h2_cc3aab(icc)
     h3 = h3_cc3aab(icc)
     !  3: (++-|++-)
     ! 12: (--+|--+)
     !tcc3(p1,p2,p3,h1,h2,h3) = +cic(idet,1)
     tcc3(p2,p1,p3,h1,h2,h3) = -tcc3(p1,p2,p3,h1,h2,h3)
     tcc3(p1,p2,p3,h2,h1,h3) = -tcc3(p1,p2,p3,h1,h2,h3)
     tcc3(p2,p1,p3,h2,h1,h3) = +tcc3(p1,p2,p3,h1,h2,h3)
  end do

end subroutine tdcc_filltcc3aab
!######################################################################
