!######################################################################
subroutine tdcc_mkden2(cic,den2)

  use mod_const,only : czero
  use mod_ormas,only : nact
  use mod_cc,only : norb1,cc_rank,den1s,den2s
  use mod_cc,only : den1_noref,den2_noref,optcc,bcc,optbcc

  implicit none
  complex(kind(0d0)),intent(in) :: cic(1:*)
  complex(kind(0d0)),intent(out) :: den2(1:*)
!#####
!DEBUG
integer :: p1,p2,q1,q2
!complex(kind(0d0)),allocatable :: den2s_pathak(:,:,:,:,:)
!DEBUG
!#####

  call tdcc_getcc(cic)
  den1s(1:nact,1:nact,1:2) = czero
  den2s(1:nact,1:nact,1:nact,1:nact,1:6) = czero

  if (cc_rank == 2 .and. optcc) then
     call ccd_den1p_main(den1s)
     call ccd_den2p_main(den2s)
  else if (cc_rank == 2 .and. bcc) then
     call bccd_den1p_main(den1s)
     call bccd_den2p_main(den2s)
  else if (cc_rank == 2 .and. optbcc) then
     call bccd_den1p_main(den1s)
     call bccd_den2p_main(den2s)
  else if (cc_rank == 2) then
     call ccsd_den1_main(den1s)
     call ccsd_den2_main(den2s)
  else if (cc_rank == 3 .and. optcc) then
     call ccdt_den1p_main(den1s)
     call ccdt_den2p_main(den2s)
  else if (cc_rank == 3 .and. bcc) then
     call bccdt_den1p_main(den1s)
     call bccdt_den2p_main(den2s)
  else if (cc_rank == 3 .and. optbcc) then
     call bccdt_den1p_main(den1s)
     call bccdt_den2p_main(den2s)
  else
     stop 'tdcc_den1 nyi.'
  end if

  den1_noref = den1s
  den2_noref = den2s

  ! tdcc_mkden2_ref2 uses den1s without reference contributions.
  ! so tdcc_mkden2_ref2 should be called before tdcc_mkden2_ref1.
  call tdcc_mkden2_ref2(den1s,den2s)

  !#####
  !OLD
  !call tdcc_mkden2_spac2(den2s,den2)
  !NEW
  call tdcc_mkden2_full(den2s,den2_noref)
  call tdcc_mkden2_aaclean(den2_noref)

!debug
!  do p1=1,nact
!  do p2=1,nact
!     do q1=1,nact
!     do q2=1,nact
!        write(6,"('tdcc_mkden2: ',4i5,6f20.15)") p1,p2,q1,q2,&
!             den2_noref(p1,p2,q1,q2,1),den2_noref(p1,p2,q1,q2,3),&
!             den2_noref(p1,p2,q1,q2,3)-den2_noref(p1,p2,q2,q1,3)
!     end do
!     end do
!  end do
!  end do
!  stop
!debug

  call tdcc_mkden2_spac2_new(den2_noref,den2)
  !#####

!  call tdcc_mkden2_ref1(den1s)
!  call tdcc_mkden2_spac1(den1s,den1)

end subroutine tdcc_mkden2
!######################################################################
subroutine tdcc_mkden2_ref2(den1s,den2s)

  use mod_bas,only : smul
  use mod_const,only : runit
  use mod_ormas,only : nact,den2_abonly
  use mod_cc,only : norb1,cc_code

  implicit none
  complex(kind(0d0)),intent(in) :: den1s(1:nact,1:nact)
  complex(kind(0d0)),intent(inout) :: den2s(1:nact,1:nact,1:nact,1:nact,1:6)

  integer :: iact,jact,kact,lact,aact,bact,cact,dact

  ! rho^{kl}_{ij}
  ! aaaa part
  if (.not.(smul==1.and.den2_abonly)) then
     do iact = 1,norb1
     do jact = 1,norb1
     do kact = 1,norb1
     do lact = 1,norb1
        if (kact == iact) den2s(kact,lact,iact,jact,1) &
                        = den2s(kact,lact,iact,jact,1) + den1s(lact,jact)
        if (lact == iact) den2s(kact,lact,iact,jact,1) &
                        = den2s(kact,lact,iact,jact,1) - den1s(kact,jact)
        if (kact == jact) den2s(kact,lact,iact,jact,1) &
                        = den2s(kact,lact,iact,jact,1) - den1s(lact,iact)
        if (lact == jact) den2s(kact,lact,iact,jact,1) &
                        = den2s(kact,lact,iact,jact,1) + den1s(kact,iact)
        if (kact == jact .and. lact == iact) den2s(kact,lact,iact,jact,1) &
                                           = den2s(kact,lact,iact,jact,1) - runit
        if (kact == iact .and. lact == jact) den2s(kact,lact,iact,jact,1) &
                                           = den2s(kact,lact,iact,jact,1) + runit
     end do
     end do
     end do
     end do
  end if
  ! abab part
  do iact = 1,norb1
  do jact = 1,norb1
  do kact = 1,norb1
  do lact = 1,norb1
     if (kact == iact) den2s(kact,lact,iact,jact,3) &
                     = den2s(kact,lact,iact,jact,3) + den1s(lact,jact)
     if (lact == jact) den2s(kact,lact,iact,jact,3) &
                     = den2s(kact,lact,iact,jact,3) + den1s(kact,iact)
     if (kact == iact .and. lact == jact) then
        den2s(kact,lact,iact,jact,3) &
      = den2s(kact,lact,iact,jact,3) + runit
     end if
  end do
  end do
  end do
  end do
  ! abba part: in fact not needed...
  do iact = 1,norb1
  do jact = 1,norb1
  do kact = 1,norb1
  do lact = 1,norb1
     if (lact == iact) den2s(kact,lact,iact,jact,5) &
                     = den2s(kact,lact,iact,jact,5) - den1s(kact,jact)
     if (kact == jact) den2s(kact,lact,iact,jact,5) &
                     = den2s(kact,lact,iact,jact,5) - den1s(lact,iact)
     if (kact == jact .and. lact == iact) then
        den2s(kact,lact,iact,jact,5) &
      = den2s(kact,lact,iact,jact,5) - runit
     end if
  end do
  end do
  end do
  end do

  ! rho^{ka}_{ij}
  ! aaaa part
  if (.not.(smul==1.and.den2_abonly)) then
     do iact = 1,norb1
     do jact = 1,norb1
     do kact = 1,norb1
     do aact = norb1 + 1,nact
        if (kact == iact) then
           den2s(kact,aact,iact,jact,1) &
         = den2s(kact,aact,iact,jact,1) + den1s(aact,jact)
        end if
        if (kact == jact) then
           den2s(kact,aact,iact,jact,1) &
         = den2s(kact,aact,iact,jact,1) - den1s(aact,iact)
        end if
     end do
     end do
     end do
     end do
  end if
  ! abab/abba part
  do iact = 1,norb1
  do jact = 1,norb1
  do kact = 1,norb1
  do aact = norb1 + 1,nact
     if (kact == iact) then
        den2s(kact,aact,iact,jact,3) &
      = den2s(kact,aact,iact,jact,3) + den1s(aact,jact)
     end if
     if (kact == jact) then
        den2s(kact,aact,iact,jact,5) &
      = den2s(kact,aact,iact,jact,5) - den1s(aact,iact)
     end if
  end do
  end do
  end do
  end do

  ! rho^{ij}_{ka}
  ! aaaa part
  if (.not.(smul==1.and.den2_abonly)) then
     do iact = 1,norb1
     do jact = 1,norb1
     do kact = 1,norb1
     do aact = norb1 + 1,nact
        if (kact == iact) then
           den2s(iact,jact,kact,aact,1) &
         = den2s(iact,jact,kact,aact,1) + den1s(jact,aact)
        end if
        if (kact == jact) then
           den2s(iact,jact,kact,aact,1) &
         = den2s(iact,jact,kact,aact,1) - den1s(iact,aact)
        end if
     end do
     end do
     end do
     end do
  end if
  ! abab/abba part
  do iact = 1,norb1
  do jact = 1,norb1
  do kact = 1,norb1
  do aact = norb1 + 1,nact
     if (kact == iact) then
        den2s(iact,jact,kact,aact,3) &
      = den2s(iact,jact,kact,aact,3) + den1s(jact,aact)
     end if
     if (kact == jact) then
        den2s(iact,jact,kact,aact,5) &
      = den2s(iact,jact,kact,aact,5) - den1s(iact,aact)
     end if
  end do
  end do
  end do
  end do

  ! rho^{jb}_{ai}
  ! aaaa part
  if (.not.(smul==1.and.den2_abonly)) then
     do iact = 1,norb1
     do jact = 1,norb1
        if (iact == jact) then
           do aact = norb1 + 1,nact
           do bact = norb1 + 1,nact
              den2s(jact,bact,aact,iact,1) &
            = den2s(jact,bact,aact,iact,1) - den1s(bact,aact)
           end do
           end do
        end if
     end do
     end do
  end if
  ! abba part
  do iact = 1,norb1
  do jact = 1,norb1
     if (iact == jact) then
        do aact = norb1 + 1,nact
        do bact = norb1 + 1,nact
           den2s(jact,bact,aact,iact,5) &
         = den2s(jact,bact,aact,iact,5) - den1s(bact,aact)
        end do
        end do
     end if
  end do
  end do

end subroutine tdcc_mkden2_ref2
!######################################################################
subroutine tdcc_mkden2_full(den2s,den2_noref)

  use mod_bas,only : smul,mval
  use mod_ormas,only : nact,den2_abonly,ncore
  use mod_cc,only : norb1,cc_code,nbiort

  implicit none
  complex(kind(0d0)),intent(in) :: den2s(1:nact,1:nact,1:nact,1:nact,1:6)
  complex(kind(0d0)),intent(inout) :: den2_noref(1:nact,1:nact,1:nact,1:nact,1:6)

  integer :: pact,qact,ract,sact
  integer :: iact,jact,kact,lact
  integer :: aact,bact,cact,dact
  integer :: p1act,p2act,q1act,q2act
  complex(kind(0d0)) :: tmp

  den2_noref = den2s

  if (.not.(smul==1.and.den2_abonly)) then
     do iact = 1,norb1
     do jact = 1,norb1
     do kact = 1,norb1
     do aact = norb1 + 1,nact
        den2_noref(aact,iact,jact,kact,1) = -den2_noref(iact,aact,jact,kact,1)
        den2_noref(iact,jact,aact,kact,1) = -den2_noref(iact,jact,kact,aact,1)
     end do
     end do
     end do
     end do
  end if
  do iact = 1,norb1
  do jact = 1,norb1
  do kact = 1,norb1
  do aact = norb1 + 1,nact
     den2_noref(aact,iact,jact,kact,3) = -den2_noref(iact,aact,jact,kact,5)
     den2_noref(aact,iact,jact,kact,5) = -den2_noref(iact,aact,jact,kact,3)
     den2_noref(iact,jact,aact,kact,3) = -den2_noref(iact,jact,kact,aact,5)
     den2_noref(iact,jact,aact,kact,5) = -den2_noref(iact,jact,kact,aact,3)
  end do
  end do
  end do
  end do

  if (.not.(smul==1.and.den2_abonly)) then
     do iact = 1,norb1
     do jact = 1,norb1
     do aact = norb1 + 1,nact
     do bact = norb1 + 1,nact
        den2_noref(aact,iact,bact,jact,1) = -den2_noref(iact,aact,bact,jact,1)
        den2_noref(aact,iact,jact,bact,1) =  den2_noref(iact,aact,bact,jact,1)
        den2_noref(iact,aact,jact,bact,1) =  den2_noref(aact,iact,bact,jact,1)
     end do
     end do
     end do
     end do
  end if
  do iact = 1,norb1
  do jact = 1,norb1
  do aact = norb1 + 1,nact
  do bact = norb1 + 1,nact
     den2_noref(aact,iact,bact,jact,3) = -den2_noref(iact,aact,bact,jact,5)
     den2_noref(aact,iact,bact,jact,5) = -den2_noref(iact,aact,bact,jact,3)
     den2_noref(aact,iact,jact,bact,3) =  den2_noref(iact,aact,bact,jact,3)
     den2_noref(aact,iact,jact,bact,5) =  den2_noref(iact,aact,bact,jact,5)
     den2_noref(iact,aact,jact,bact,3) =  den2_noref(aact,iact,bact,jact,3)
     den2_noref(iact,aact,jact,bact,5) =  den2_noref(aact,iact,bact,jact,5)
  end do
  end do
  end do
  end do

  if (.not.(smul==1.and.den2_abonly)) then
     do iact = 1,norb1
     do aact = norb1 + 1,nact
     do bact = norb1 + 1,nact
     do cact = norb1 + 1,nact
        den2_noref(iact,aact,bact,cact,1) = -den2_noref(aact,iact,bact,cact,1)
        den2_noref(aact,bact,iact,cact,1) = -den2_noref(aact,bact,cact,iact,1)
     end do
     end do
     end do
     end do
  end if
  do iact = 1,norb1
  do aact = norb1 + 1,nact
  do bact = norb1 + 1,nact
  do cact = norb1 + 1,nact
     den2_noref(iact,aact,bact,cact,3) = -den2_noref(aact,iact,bact,cact,5)
     den2_noref(iact,aact,bact,cact,5) = -den2_noref(aact,iact,bact,cact,3)
     den2_noref(aact,bact,iact,cact,3) = -den2_noref(aact,bact,cact,iact,5)
     den2_noref(aact,bact,iact,cact,5) = -den2_noref(aact,bact,cact,iact,3)
  end do
  end do
  end do
  end do

  if (smul==1.and.den2_abonly) then
     den2_noref(:,:,:,:,1) = 0d0
     do p1act = 1,nact
     do p2act = 1,nact
     do q1act = 1,nact
     do q2act = 1,nact
        if (mval(ncore+p1act)+mval(ncore+p2act).ne.mval(ncore+q1act)+mval(ncore+q2act) .or. &
             p1act==p2act .or. q1act==q2act) cycle
        den2_noref(p1act,p2act,q1act,q2act,1) = &
        den2_noref(p1act,p2act,q1act,q2act,3) - &
        den2_noref(p1act,p2act,q2act,q1act,3)
     end do
     end do
     end do
     end do
  end if

  den2_noref(:,:,:,:,2) = den2_noref(:,:,:,:,1)
  den2_noref(:,:,:,:,4) = den2_noref(:,:,:,:,3)
  den2_noref(:,:,:,:,6) = den2_noref(:,:,:,:,5)

end subroutine tdcc_mkden2_full
!######################################################################
subroutine tdcc_mkden2_spac2(den2s,den2)

  use mod_const,only : ctwo,chalf
  use mod_ormas,only : nact
  use mod_cc,only : norb1,nbiort

  implicit none
  complex(kind(0d0)),intent(inout) :: den2s(1:nact,1:nact,1:nact,1:nact,1:6)
  complex(kind(0d0)),intent(out) :: den2(1:nact,1:nact,1:nact,1:nact)

  integer :: pact,qact,ract,sact
  integer :: iact,jact,kact,lact
  integer :: aact,bact,cact,dact
  complex(kind(0d0)) :: tmp
  complex(kind(0d0)),allocatable :: den2t(:,:,:,:,:)
  complex(kind(0d0)),allocatable :: den2u(:,:,:,:,:)

  allocate(den2t(1:nact,1:nact,1:nact,1:nact,1:6))
  allocate(den2u(1:nact,1:nact,1:nact,1:nact,1:6))
  den2t = den2s

  ! den2  = 1/2 (den2t + den2t^+)
  ! den2t^{rs}_{pq} = 2 den2s^{rs}_{pq}(aaaa) + den2s^{rs}_{pq}(aabb) + den2s^{sr}_{qp}(aabb)
  
  do iact = 1,norb1
  do jact = 1,norb1
  do kact = 1,norb1
  do lact = 1,norb1
     den2t(iact,jact,kact,lact,1) = (den2t(iact,jact,kact,lact,1) &
                                      +  den2t(iact,jact,kact,lact,3)) * ctwo

     den2u(iact,jact,kact,lact,1) = den2s(iact,jact,kact,lact,1)
     den2u(iact,jact,kact,lact,3) = den2s(iact,jact,kact,lact,3)
     den2u(iact,jact,kact,lact,5) = den2s(iact,jact,kact,lact,5)
  end do
  end do
  end do
  end do

  do iact = 1,norb1
  do jact = 1,norb1
  do kact = 1,norb1
  do aact = norb1 + 1,nact
     tmp = den2t(iact,aact,jact,kact,1)
     den2t(iact,aact,jact,kact,1) = (tmp + den2t(iact,aact,jact,kact,3)) * (+ ctwo)
     den2t(aact,iact,jact,kact,1) = (tmp + den2t(iact,aact,jact,kact,5)) * (- ctwo)

     den2u(iact,aact,jact,kact,1) = +den2s(iact,aact,jact,kact,1)
     den2u(iact,aact,jact,kact,3) = +den2s(iact,aact,jact,kact,3)
     den2u(iact,aact,jact,kact,5) = +den2s(iact,aact,jact,kact,5)
     den2u(aact,iact,jact,kact,1) = -den2s(iact,aact,jact,kact,1)
     den2u(aact,iact,jact,kact,3) = -den2s(iact,aact,jact,kact,5)
     den2u(aact,iact,jact,kact,5) = -den2s(iact,aact,jact,kact,3)

     tmp = den2t(iact,jact,kact,aact,1)
     den2t(iact,jact,kact,aact,1) = (tmp + den2t(iact,jact,kact,aact,3)) * (+ ctwo)
     den2t(iact,jact,aact,kact,1) = (tmp + den2t(iact,jact,kact,aact,5)) * (- ctwo)

     den2u(iact,jact,kact,aact,1) = +den2s(iact,jact,kact,aact,1)
     den2u(iact,jact,kact,aact,3) = +den2s(iact,jact,kact,aact,3)
     den2u(iact,jact,kact,aact,5) = +den2s(iact,jact,kact,aact,5)
     den2u(iact,jact,aact,kact,1) = -den2s(iact,jact,kact,aact,1)
     den2u(iact,jact,aact,kact,3) = -den2s(iact,jact,kact,aact,5)
     den2u(iact,jact,aact,kact,5) = -den2s(iact,jact,kact,aact,3)
  end do
  end do
  end do
  end do

  do iact = 1,norb1
  do jact = 1,norb1
  do aact = norb1 + 1,nact
  do bact = norb1 + 1,nact
     den2t(iact,jact,aact,bact,1) = (den2t(iact,jact,aact,bact,1) &
                                      +  den2t(iact,jact,aact,bact,3)) * ctwo
     den2t(aact,bact,iact,jact,1) = (den2t(aact,bact,iact,jact,1) &
                                      +  den2t(aact,bact,iact,jact,3)) * ctwo

     den2u(iact,jact,aact,bact,1) = den2s(iact,jact,aact,bact,1)
     den2u(iact,jact,aact,bact,3) = den2s(iact,jact,aact,bact,3)
     den2u(iact,jact,aact,bact,5) = den2s(iact,jact,aact,bact,5)
     den2u(aact,bact,iact,jact,1) = den2s(aact,bact,iact,jact,1)
     den2u(aact,bact,iact,jact,3) = den2s(aact,bact,iact,jact,3)
     den2u(aact,bact,iact,jact,5) = den2s(aact,bact,iact,jact,5)

     tmp = den2t(iact,aact,bact,jact,1)
     den2t(iact,aact,bact,jact,1) = (tmp + den2t(iact,aact,bact,jact,3)) * (+ ctwo)
     den2t(aact,iact,bact,jact,1) = (tmp + den2t(iact,aact,bact,jact,5)) * (- ctwo)

     den2u(iact,aact,bact,jact,1) = +den2s(iact,aact,bact,jact,1)
     den2u(iact,aact,bact,jact,3) = +den2s(iact,aact,bact,jact,3)
     den2u(iact,aact,bact,jact,5) = +den2s(iact,aact,bact,jact,5)
     den2u(aact,iact,bact,jact,1) = -den2s(iact,aact,bact,jact,1)
     den2u(aact,iact,bact,jact,3) = -den2s(iact,aact,bact,jact,5)
     den2u(aact,iact,bact,jact,5) = -den2s(iact,aact,bact,jact,3)

     den2t(aact,iact,jact,bact,1) =  den2t(iact,aact,bact,jact,1)
     den2t(iact,aact,jact,bact,1) =  den2t(aact,iact,bact,jact,1)

     den2u(aact,iact,jact,bact,1) =  den2u(iact,aact,bact,jact,1)
     den2u(aact,iact,jact,bact,3) =  den2u(iact,aact,bact,jact,3)
     den2u(aact,iact,jact,bact,5) =  den2u(iact,aact,bact,jact,5)
     den2u(iact,aact,jact,bact,1) =  den2u(aact,iact,bact,jact,1)
     den2u(iact,aact,jact,bact,3) =  den2u(aact,iact,bact,jact,3)
     den2u(iact,aact,jact,bact,5) =  den2u(aact,iact,bact,jact,5)
  end do
  end do
  end do
  end do

  do iact = 1,norb1
  do aact = norb1 + 1,nact
  do bact = norb1 + 1,nact
  do cact = norb1 + 1,nact
     tmp = den2t(aact,iact,bact,cact,1)
     den2t(aact,iact,bact,cact,1) = (tmp + den2t(aact,iact,bact,cact,3)) * (+ ctwo)
     den2t(iact,aact,bact,cact,1) = (tmp + den2t(aact,iact,bact,cact,5)) * (- ctwo)

     den2u(aact,iact,bact,cact,1) = +den2s(aact,iact,bact,cact,1)
     den2u(aact,iact,bact,cact,3) = +den2s(aact,iact,bact,cact,3)
     den2u(aact,iact,bact,cact,5) = +den2s(aact,iact,bact,cact,5)
     den2u(iact,aact,bact,cact,1) = -den2s(aact,iact,bact,cact,1)
     den2u(iact,aact,bact,cact,3) = -den2s(aact,iact,bact,cact,5)
     den2u(iact,aact,bact,cact,5) = -den2s(aact,iact,bact,cact,3)

     tmp = den2t(aact,bact,cact,iact,1)
     den2t(aact,bact,cact,iact,1) = (tmp + den2t(aact,bact,cact,iact,3)) * (+ ctwo)
     den2t(aact,bact,iact,cact,1) = (tmp + den2t(aact,bact,cact,iact,5)) * (- ctwo)

     den2u(aact,bact,cact,iact,1) = +den2s(aact,bact,cact,iact,1)
     den2u(aact,bact,cact,iact,3) = +den2s(aact,bact,cact,iact,3)
     den2u(aact,bact,cact,iact,5) = +den2s(aact,bact,cact,iact,5)
     den2u(aact,bact,iact,cact,1) = -den2s(aact,bact,cact,iact,1)
     den2u(aact,bact,iact,cact,3) = -den2s(aact,bact,cact,iact,5)
     den2u(aact,bact,iact,cact,5) = -den2s(aact,bact,cact,iact,3)
  end do
  end do
  end do
  end do

  do aact = norb1 + 1,nact
  do bact = norb1 + 1,nact
  do cact = norb1 + 1,nact
  do dact = norb1 + 1,nact
     den2t(aact,bact,cact,dact,1) = (den2t(aact,bact,cact,dact,1) &
                                      +  den2t(aact,bact,cact,dact,3)) * ctwo

     den2u(aact,bact,cact,dact,1) = den2s(aact,bact,cact,dact,1)
     den2u(aact,bact,cact,dact,3) = den2s(aact,bact,cact,dact,3)
     den2u(aact,bact,cact,dact,5) = den2s(aact,bact,cact,dact,5)
  end do
  end do
  end do
  end do

  if (nbiort == 1) then
     do pact = 1,nact
     do qact = 1,nact
     do ract = 1,nact
     do sact = 1,nact
        tmp =       den2t(pact,ract,qact,sact,1) &
            + conjg(den2t(qact,sact,pact,ract,1))
        den2(pact,qact,ract,sact) = tmp * chalf
     end do
     end do
     end do
     end do
  else
     do pact = 1,nact
     do qact = 1,nact
     do ract = 1,nact
     do sact = 1,nact
        den2(pact,qact,ract,sact) = den2t(pact,ract,qact,sact,1)
     end do
     end do
     end do
     end do
  end if

!  !DEBUG
!  write(6,"('den2cc-R')")
! ! do pact = 1,norb1
! ! do qact = 1,norb1
! !    do ract = 1,norb1
! !    do sact = 1,norb1
!  do pact = 1,nact
!  do qact = 1,nact
!  do ract = 1,nact
!  do sact = 1,nact
!     write(6,"(4i5,e20.10)") pact,qact,ract,sact,dble(den2(pact,qact,ract,sact))
!  end do
!  end do
!  end do
!  end do
!  stop
!  !DEBUG

! den2s = den2t * chalf
  den2s = den2u
  den2s(:,:,:,:,2) = den2s(:,:,:,:,1)
  den2s(:,:,:,:,4) = den2s(:,:,:,:,3)
  den2s(:,:,:,:,6) = den2s(:,:,:,:,5)

  deallocate(den2u)
  deallocate(den2t)

end subroutine tdcc_mkden2_spac2
!######################################################################
subroutine tdcc_mkden2_spac2_new(den2s,den2)

  use mod_const,only : ctwo,chalf
  use mod_bas,only : mval,smul
  use mod_ormas,only : nact,ncore,nelact
  use mod_cc,only : norb1,nbiort

  implicit none
  complex(kind(0d0)),intent(in) :: den2s(1:nact,1:nact,1:nact,1:nact,1:6)
  complex(kind(0d0)),intent(out) :: den2(1:nact,1:nact,1:nact,1:nact)

  integer :: pact,qact,ract,sact
  integer :: iact,jact,kact,lact
  integer :: aact,bact,cact,dact
  complex(kind(0d0)) :: tmp

  ! den2  = 1/2 (den2t + den2t^+)
  ! den2t^{rs}_{pq} = 2 den2s^{rs}_{pq}(aaaa) + den2s^{rs}_{pq}(aabb) + den2s^{sr}_{qp}(aabb)
  
  if (nbiort == 1) then
     if (smul==1.and.nelact(1)==nelact(2)) then
        do pact = 1,nact
        do qact = 1,nact
        do ract = 1,nact
        do sact = 1,nact
           if (mval(ncore+pact)+mval(ncore+ract) &
          .ne. mval(ncore+qact)+mval(ncore+sact)) cycle
           den2(pact,qact,ract,sact) = chalf*( &
                + den2s(pact,ract,qact,sact,1) + den2s(pact,ract,qact,sact,1) &
                + den2s(pact,ract,qact,sact,3) + den2s(ract,pact,sact,qact,3) &
          + conjg(den2s(qact,sact,pact,ract,1) + den2s(qact,sact,pact,ract,1) &
                + den2s(qact,sact,pact,ract,3) + den2s(sact,qact,ract,pact,3)))
        end do
        end do
        end do
        end do
     else if (nelact(2)==0) then
        do pact = 1,nact
        do qact = 1,nact
        do ract = 1,nact
        do sact = 1,nact
           if (mval(ncore+pact)+mval(ncore+ract) &
          .ne. mval(ncore+qact)+mval(ncore+sact)) cycle
           den2(pact,qact,ract,sact) = chalf*( &
                + den2s(pact,ract,qact,sact,1) &
          + conjg(den2s(qact,sact,pact,ract,1)))
        end do
        end do
        end do
        end do
     else
        stop 'tdcc_mkden2_spac2_new: A=B.or.B=0 only'
     end if
  else
     do pact = 1,nact
     do qact = 1,nact
     do ract = 1,nact
     do sact = 1,nact
        if (mval(ncore+pact)+mval(ncore+ract) &
       .ne. mval(ncore+qact)+mval(ncore+sact)) cycle
        den2(pact,qact,ract,sact) = &
             + den2s(pact,ract,qact,sact,1) + den2s(pact,ract,qact,sact,1) &
             + den2s(pact,ract,qact,sact,3) + den2s(ract,pact,sact,qact,4)
     end do
     end do
     end do
     end do
  end if

!  !DEBUG
!  write(6,"('den2cc-R')")
! ! do pact = 1,norb1
! ! do qact = 1,norb1
! !    do ract = 1,norb1
! !    do sact = 1,norb1
!  do pact = 1,nact
!  do qact = 1,nact
!  do ract = 1,nact
!  do sact = 1,nact
!     write(6,"(4i5,e20.10)") pact,qact,ract,sact,dble(den2(pact,qact,ract,sact))
!  end do
!  end do
!  end do
!  end do
!  stop
!  !DEBUG

end subroutine tdcc_mkden2_spac2_new
!######################################################################
subroutine tdcc_mkden2_aaclean(den2)

  !antisymmetrize aaaa part of 2RDM

  use mod_ormas,only : nact

  implicit none
  complex(kind(0d0)),intent(inout) :: den2(1:nact,1:nact,1:nact,1:nact)

  integer :: iact,jact,kact,lact

  do iact = 1,nact
     do kact = 1,nact
        do lact = 1,nact
           den2(iact,iact,kact,lact) = 0d0
        end do
     end do
  end do
  do iact = 1,nact
     do jact = 1,nact
        do kact = 1,nact
           den2(iact,jact,kact,kact) = 0d0
        end do
     end do
  end do

  do iact = 1,nact
     do jact = 1,iact-1
        do kact = 1,nact
           do lact = 1,kact-1
             !den2(iact,jact,kact,lact) = +den2(iact,jact,kact,lact)
              den2(iact,jact,lact,kact) = -den2(iact,jact,kact,lact)
              den2(jact,iact,kact,lact) = -den2(iact,jact,kact,lact)
              den2(jact,iact,lact,kact) = +den2(iact,jact,kact,lact)
           end do
        end do
     end do
  end do

end subroutine tdcc_mkden2_aaclean
!######################################################################
