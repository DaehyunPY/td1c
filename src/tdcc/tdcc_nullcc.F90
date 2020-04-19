!######################################################################
subroutine tdcc_nullcc(cic)

  use, intrinsic :: iso_c_binding
  use mod_cc, only : cc_rank
  implicit none
  complex(c_double_complex), intent(out) :: cic(1:*)

  call tdcc_nulltcc0(cic)
  if (cc_rank >= 1) call tdcc_nulltcc1(cic)
  if (cc_rank >= 2) call tdcc_nulltcc2(cic)
  if (cc_rank >= 3) call tdcc_nulltcc3(cic)

  call tdcc_nullgcc0(g0out,cic)
  if (cc_rank >= 1) call tdcc_nullgcc1(cic)
  if (cc_rank >= 2) call tdcc_nullgcc2(cic)
  if (cc_rank >= 3) call tdcc_nullgcc3(cic)

end subroutine tdcc_nullcc
!######################################################################
subroutine tdcc_nulltcc0(cic)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : ndetx
  use mod_cc, only : map_cc0
  implicit none
  complex(c_double_complex), intent(out) :: cic(1:ndetx,1:2)

  cic(map_cc0(1,3),1) = tcc0

end subroutine tdcc_nulltcc0
!######################################################################
subroutine tdcc_nullgcc0(cic)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : ndetx
  use mod_cc, only : map_cc0
  implicit none
  complex(c_double_complex), intent(out) :: cic(1:ndetx,1:2)

  cic(map_cc0(1,3),2) = gcc0

end subroutine tdcc_nullgcc0
!######################################################################
subroutine tdcc_nulltcc1(cic)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nact,ndetx
  use mod_cc, only : norb1,ncc1a
  implicit none
  complex(c_double_complex), intent(out) :: cic(1:ndetx,1:2)
  integer(c_int) :: icc,idet

  do icc = 1, ncc1a
     idet = map_cc1a(icc,3)
     cic(idet,1) = 0d0
  end do

end subroutine tdcc_nulltcc1
!######################################################################
subroutine tdcc_nullgcc1(cic)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nact,ndetx
  use mod_cc, only : norb1,ncc1a
  implicit none
  complex(c_double_complex), intent(out) :: cic(1:ndetx,1:2)
  integer(c_int) :: icc,idet

  do icc = 1, ncc1a
     idet = map_cc1a(icc,3)
     cic(idet,2) = 0d0
  end do

end subroutine tdcc_nullgcc1
!######################################################################
subroutine tdcc_nulltcc2(cic)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nact,ndetx
  use mod_cc, only : norb1,ncc2aa,ncc2ab,map_cc2aa,map_cc2ab
  implicit none
  complex(c_double_complex), intent(out) :: cic(1:ndetx,1:2)
  integer(c_int) :: icc,idet

  ! (p+p+|h+h+)
  do icc = 1, ncc2aa
     idet = map_cc2aa(icc,3)
     cic(idet,1) = 0d0
  end do

  ! (p+p-|h+h-)
  do icc = 1, ncc2ab
     idet = map_cc2ab(icc,3)
     cic(idet,1) = 0d0
  end do

end subroutine tdcc_nulltcc2
!######################################################################
subroutine tdcc_nullgcc2(cic)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nact,ndetx
  use mod_cc, only : norb1,ncc2aa,ncc2ab,map_cc2aa,map_cc2ab
  implicit none
  complex(c_double_complex), intent(out) :: cic(1:ndetx,1:2)
  integer(c_int) :: icc,idet

  ! (h+h+|p+p+)
  do icc = 1, ncc2aa
     idet = map_cc2aa(icc,3)
     cic(idet,2) = 0d0
  end do

  ! (h+h-|p+p-)
  do icc = 1, ncc2ab
     idet = map_cc2ab(icc,3)
     cic(idet,2) = 0d0
  end do

end subroutine tdcc_nullgcc2
!######################################################################
subroutine tdcc_nulltcc3(cic)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nact,ndetx
  use mod_cc, only : norb1,ncc3aaa,ncc3aab,map_cc3aaa,map_cc3aab
  implicit none
  complex(c_double_complex), intent(out) :: cic(1:ndetx,1:2)
  integer(c_int) :: icc,idet

  ! (+++|+++)
  do icc = 1, ncc3aaa
     idet = map_cc3aaa(icc,3)
     cic(idet,1) = 0d0
  end do

  ! (++-|++-)
  do icc = 1, ncc3aab
     idet = map_cc3aab(icc,3)
     cic(idet,1) = 0d0
  end do

end subroutine tdcc_nulltcc3
!######################################################################
subroutine tdcc_nullgcc3(cic)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nact,ndetx
  use mod_cc, only : norb1,ncc3aaa,ncc3aab,map_cc3aaa,map_cc3aab
  implicit none
  complex(c_double_complex), intent(out) :: cic(1:ndetx,1:2)
  integer(c_int) :: icc,idet

  ! (+++|+++)
  do icc = 1, ncc3aaa
     idet = map_cc3aaa(icc,3)
     cic(idet,2) = 0d0
  end do

  ! (++-|++-)
  do icc = 1, ncc3aab
     idet = map_cc3aab(icc,3)
     cic(idet,2) = 0d0
  end do

end subroutine tdcc_nullgcc3
!######################################################################
