!################################################################################
subroutine tdcc_ccsolve(int1e,int2e,cic)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : ndetx
  use mod_cc, only : cc_rank,dot1,tonly
  use mod_cc, only : thrtamp,thrgamp,cc_maxcyc
  use mod_cc, only : fock,int2x

  implicit none
  !--------------------------------------------------------------------
  complex(c_double_complex), intent(in) :: int1e(*),int2e(*)
  complex(c_double_complex), intent(out) :: cic(ndetx,2)
  !--------------------------------------------------------------------
  integer(c_long) :: icyc
  real(c_double) :: ene,ene1,dt1,dt2(2),dt3(2)
  complex(c_double_complex), allocatable :: hcic(:,:)
  real(c_double), external :: tdcc_enec
  complex(c_double_complex), external :: tdcc_ene2

  allocate(hcic(1:ndetx,1:2))

  dt1 = 1d0
  dt2 = 1d0
  dt3 = 1d0
  ene1 = 1d0

  call tdcc_mkint1x(int1e, int2e, fock)
  call tdcc_mkint2x(int2e, int2x)

  !##### INITIALIZE T amplitudes #####
  cic(1:ndetx,1) = 0d0

  do icyc = 1, cc_maxcyc
     hcic = 0d0
     call tdcc_hcc12(int1e,int2e,cic,hcic)

     ene = tdcc_enec(int1e,fock)+dble(tdcc_ene2())
     write(6,"('tdcc_solve_t: ',i5,f20.10,6e20.10)") icyc,ene,ene-ene1,dt1,dt2(1:2),dt3(1:2)
     if ((cc_rank < 1 .or. .not.dot1 .or. abs(dt1)< thrtamp) .and. &
         (cc_rank < 2 .or. (abs(dt2(1)) < thrtamp .and. abs(dt2(2)) < thrtamp)) .and. &
         (cc_rank < 3 .or. (abs(dt3(1)) < thrtamp .and. abs(dt3(2)) < thrtamp))) exit
     ene1 = ene

     if (dot1) then
        call tdcc_solve1(cic,hcic,dt1)
        hcic = 0d0
        call tdcc_hcc12(int1e,int2e,cic,hcic)
     end if

     if (cc_rank >=2) call tdcc_solve2(cic,hcic,dt2)

     if (cc_rank >=3) then
        hcic = 0d0
        call tdcc_hcc12(int1e,int2e,cic,hcic)
        call tdcc_solve3(cic,hcic,dt3)
     end if
  end do

  dt1 = 1d0
  dt2 = 1d0
  dt3 = 1d0
  if (.not. tonly) then
     !##### INITIALIZE L amplitudes #####
     cic(1:ndetx,2) = 0d0

     do icyc = 1, cc_maxcyc
        write(6,"('tdcc_solve_l: ',i5,40x,5e20.10)") icyc,dt1,dt2(1:2),dt3(1:2)
        if ((cc_rank < 1 .or. .not.dot1 .or. abs(dt1)< thrgamp) .and. &
            (cc_rank < 2 .or. (abs(dt2(1)) < thrgamp .and. abs(dt2(2)) < thrgamp)) .and. &
            (cc_rank < 3 .or. (abs(dt3(1)) < thrgamp .and. abs(dt3(2)) < thrgamp))) exit
     
        hcic = 0d0
        call tdcc_hcc12(int1e,int2e,cic,hcic)

        if (dot1) then
           call tdcc_solve1(cic(:,2),hcic(:,2),dt1)
           hcic = 0d0
           call tdcc_hcc12(int1e,int2e,cic,hcic)
        end if

        if (cc_rank >=2) call tdcc_solve2(cic(:,2),hcic(:,2),dt2)

        if (cc_rank >=3) then
           hcic = 0d0
           call tdcc_hcc12(int1e,int2e,cic,hcic)
           call tdcc_solve3(cic(:,2),hcic(:,2),dt3)
        end if
     end do
  end if

  deallocate(hcic)

  !debug
  !write(6,"('tdcc_ccsolve: converged amplitudes')")
  !call tdcc_print(cic)
  !stop 'for debug @ tdcc_ccsolve.'
  !debug

end subroutine tdcc_ccsolve
!######################################################################
subroutine tdcc_solve1(cic,hcic,dt1)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : ndetx
  use mod_cc, only : fock,norb1,norb2,ncc1a
  use mod_cc, only : h1_cc1a,p1_cc1a,map_cc1a

  implicit none
  complex(c_double_complex), intent(inout) :: cic(1:ndetx)
  complex(c_double_complex), intent(in) :: hcic(1:ndetx)
  real(c_double), intent(out) :: dt1
  integer(c_long) :: h1,p1,icc,idet
  complex(c_double_complex) :: tcc,rhs

  dt1=0d0
  do icc = 1, ncc1a
     h1 = h1_cc1a(icc)
     p1 = p1_cc1a(icc)
     idet = map_cc1a(icc,3)
     tcc = cic(idet)
     rhs = hcic(idet)/(fock(h1,h1,1)-fock(p1,p1,1))
     cic(idet) = tcc + rhs
     dt1 = dt1 + rhs**2
  end do
  dt1 = sqrt(dt1)/(norb1*norb2)

end subroutine tdcc_solve1
!######################################################################
subroutine tdcc_solve2(cic,hcic,dt2)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : ndetx
  use mod_cc, only : fock,norb1,norb2,ncc2aa,ncc2ab
  use mod_cc, only : h1_cc2aa,h2_cc2aa,p1_cc2aa,p2_cc2aa,map_cc2aa
  use mod_cc, only : h1_cc2ab,h2_cc2ab,p1_cc2ab,p2_cc2ab,map_cc2ab

  implicit none
  complex(c_double_complex), intent(inout) :: cic(ndetx)
  complex(c_double_complex), intent(in) :: hcic(ndetx)
  real(c_double), intent(out) :: dt2(2)
  integer(c_long) :: h1,h2,p1,p2,icc,idet
  complex(c_double_complex) :: tcc,rhs

  dt2(1)=0d0
  do icc = 1, ncc2aa
     h1 = h1_cc2aa(icc)
     h2 = h2_cc2aa(icc)
     p1 = p1_cc2aa(icc)
     p2 = p2_cc2aa(icc)
     idet = map_cc2aa(icc,3)
     tcc = cic(idet)
     rhs = hcic(idet)/(fock(h1,h1,1)+fock(h2,h2,1)-fock(p1,p1,1)-fock(p2,p2,1))
     cic(idet) = tcc + rhs
     dt2(1) = dt2(1) + rhs**2
  end do
  dt2(1) = sqrt(dt2(1))/((norb1*norb2)**2)

  dt2(2)=0d0
  do icc = 1, ncc2ab
     h1 = h1_cc2ab(icc)
     h2 = h2_cc2ab(icc)
     p1 = p1_cc2ab(icc)
     p2 = p2_cc2ab(icc)
     idet = map_cc2ab(icc,3)
     tcc = cic(idet)
     rhs = hcic(idet)/(fock(h1,h1,1)+fock(h2,h2,1)-fock(p1,p1,1)-fock(p2,p2,1))
     cic(idet) = tcc + rhs
     dt2(2) = dt2(2) + rhs**2
  end do
  dt2(2) = sqrt(dt2(2))/((norb1*norb2)**2)

end subroutine tdcc_solve2
!######################################################################
subroutine tdcc_solve3(cic,hcic,dt3)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : ndetx
  use mod_cc, only : fock,norb1,norb2,ncc3aaa,ncc3aab
  use mod_cc, only : h1_cc3aaa,h2_cc3aaa,h3_cc3aaa,p1_cc3aaa,p2_cc3aaa,p3_cc3aaa,map_cc3aaa
  use mod_cc, only : h1_cc3aab,h2_cc3aab,h3_cc3aab,p1_cc3aab,p2_cc3aab,p3_cc3aab,map_cc3aab

  implicit none
  complex(c_double_complex), intent(inout) :: cic(ndetx)
  complex(c_double_complex), intent(in) :: hcic(ndetx)
  real(c_double), intent(out) :: dt3(2)
  integer(c_long) :: h1,h2,h3,p1,p2,p3,icc,idet
  complex(c_double_complex) :: tcc,rhs

  dt3(1)=0d0
  do icc = 1, ncc3aaa
     h1 = h1_cc3aaa(icc)
     h2 = h2_cc3aaa(icc)
     h3 = h3_cc3aaa(icc)
     p1 = p1_cc3aaa(icc)
     p2 = p2_cc3aaa(icc)
     p3 = p3_cc3aaa(icc)
     idet = map_cc3aaa(icc,3)
     tcc = cic(idet)
     rhs = hcic(idet)/(fock(h1,h1,1)+fock(h2,h2,1)+fock(h3,h3,1)-fock(p1,p1,1)-fock(p2,p2,1)-fock(p3,p3,1))
     cic(idet) = tcc + rhs
     dt3(1) = dt3(1) + rhs**2
  end do
  dt3(1) = sqrt(dt3(1))/((norb1*norb2)**3)

  dt3(2)=0d0
  do icc = 1, ncc3aab
     h1 = h1_cc3aab(icc)
     h2 = h2_cc3aab(icc)
     h3 = h3_cc3aab(icc)
     p1 = p1_cc3aab(icc)
     p2 = p2_cc3aab(icc)
     p3 = p3_cc3aab(icc)
     idet = map_cc3aab(icc,3)
     tcc = cic(idet)
     rhs = hcic(idet)/(fock(h1,h1,1)+fock(h2,h2,1)+fock(h3,h3,2)-fock(p1,p1,1)-fock(p2,p2,1)-fock(p3,p3,2))
     cic(idet) = tcc + rhs
     dt3(2) = dt3(2) + rhs**2
  end do
  dt3(2) = sqrt(dt3(2))/((norb1*norb2)**3)

end subroutine tdcc_solve3
!######################################################################
!nyisubroutine tdcc_l1pt1(tcc1,gcc1)
!nyi
!nyi  use, intrinsic :: iso_c_binding
!nyi  use omp_mod
!nyi  use cc_mod, only : fock
!nyi  use wfn_mod, only : nact,norb1,norb2
!nyi
!nyi  implicit none
!nyi  complex(c_double_complex), intent(in) :: tcc1((norb1+1):nact,1:norb1,1:*)
!nyi  complex(c_double_complex), intent(out) :: gcc1(1:norb1,(norb1+1):nact,1:*)
!nyi
!nyi  integer(c_long) :: h1,p2
!nyi  do h1 = 1, norb1
!nyi  do p2 = norb1+1, nact
!nyi     gcc1(h1,p2,1) = tcc1(p2,h1,1)
!nyi  end do
!nyi  end do
!nyi
!nyiend subroutine tdcc_l1pt1
!nyi!######################################################################
!nyisubroutine tdcc_l2pt1(tcc2,gcc2)
!nyi
!nyi  use, intrinsic :: iso_c_binding
!nyi  use omp_mod
!nyi  use cc_mod, only : fock
!nyi  use wfn_mod, only : nact,norb1,norb2
!nyi
!nyi  implicit none
!nyi  complex(c_double_complex), intent(in) :: tcc2((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:*)
!nyi  complex(c_double_complex), intent(out) :: gcc2(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,1:*)
!nyi
!nyi  integer(c_long) :: h1,h2,p3,p4
!nyi  do h1 = 1, norb1
!nyi  do h2 = 1, norb1
!nyi  do p3 = norb1 + 1, nact
!nyi  do p4 = norb1 + 1, nact
!nyi     gcc2(h1,h2,p3,p4,1) = tcc2(p3,p4,h1,h2,1)
!nyi     gcc2(h1,h2,p3,p4,2) = tcc2(p3,p4,h1,h2,2)
!nyi  end do
!nyi  end do
!nyi  end do
!nyi  end do
!nyi
!nyiend subroutine tdcc_l2pt1
!nyi!######################################################################
!nyisubroutine tdcc_l3pt1(tcc3,gcc3)
!nyi
!nyi  use, intrinsic :: iso_c_binding
!nyi  use omp_mod
!nyi  use cc_mod, only : fock
!nyi  use wfn_mod, only : nact,norb1,norb2
!nyi
!nyi  implicit none
!nyi  complex(c_double_complex), intent(in) :: tcc3((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1,1:*)
!nyi  complex(c_double_complex), intent(out) :: gcc3(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:*)
!nyi
!nyi  integer(c_long) :: h1,h2,h3,p4,p5,p6
!nyi  do h1 = 1, norb1
!nyi  do h2 = 1, norb1
!nyi  do h3 = 1, norb1
!nyi  do p4 = norb1 + 1, nact
!nyi  do p5 = norb1 + 1, nact
!nyi  do p6 = norb1 + 1, nact
!nyi     gcc3(h1,h2,h3,p4,p5,p6,1) = tcc3(p4,p5,p6,h1,h2,h3,1)
!nyi     gcc3(h1,h2,h3,p4,p5,p6,2) = tcc3(p4,p5,p6,h1,h2,h3,2)
!nyi  end do
!nyi  end do
!nyi  end do
!nyi  end do
!nyi  end do
!nyi  end do
!nyi
!nyiend subroutine tdcc_l3pt1
!######################################################################
