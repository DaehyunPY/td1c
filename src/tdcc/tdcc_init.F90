!################################################################################
subroutine tdcc_init(mtot)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nelact,nact,ndetx,nsub,norb_sub,max_sub
  use mod_ormas, only : act1_ll,act1_ul,nrotaa,rotaa_mapb
  use mod_cc, only : cc_type,dot1,optcc,bcc,optbcc
  use mod_cc, only : norb1,norb2,cc_rank,fock,int2x,den1s,den2s,den1_noref,den2_noref
  use mod_cc, only : t1inp,g1inp,t2inp,g2inp,t3inp,g3inp,dt1inp,dg1inp,dt2inp,dg2inp,dt3inp,dg3inp
  use mod_cc, only : t1out,g1out,t2out,g2out,t3out,g3out,dt1out,dg1out,dt2out,dg2out,dt3out,dg3out
  use mod_cc, only : nXai,nXij,nXab,aiX,ijX,abX
  use mod_cc2
  
  implicit none
  integer(c_int), intent(in) :: mtot
  integer(c_int),external :: tdcc_spin_int2x,tdcc_spin_fock
  integer(c_int),external :: tdcc_spin_t2inp,tdcc_spin_t3inp
  integer(c_int),external :: tdcc_spin_g2inp,tdcc_spin_g3inp
  integer(c_int) :: p1,p2,irot

  if (ndetx < 2) stop "tdcc_init: ndetx < 2."
  if (nact == 0) stop "tdcc_init: nact = 0."
!  if (nsub .ne. 2) stop "tdcc_init: nsub .ne. 2."
  if (nelact(3) < 2) stop "tdcc_init: nelact < 2."
  if (nelact(1).ne.nelact(2).and.nelact(2)>0) stop "tdcc_init: nelact(1).ne.nelact(2)>0."

  norb1 = norb_sub(1)
  norb2 = norb_sub(2)
  cc_rank = max_sub(2)
  optcc = .false.
  bcc = .false.
  optbcc = .false.
  dot1 = .false.
  if (cc_type == 0) then
     optcc = .true.
  else if (cc_type == 1) then
     bcc = .true.
  else if (cc_type == 2) then
     optbcc = .true.
  else
     dot1 = .true.
  end if
  call tdcc_mkmap

  allocate(fock(1:nact,1:nact,1:2))
  allocate(int2x(1:nact,1:nact,1:nact,1:nact,1:6))
  allocate(den1s(1:nact,1:nact,1:2))
  allocate(den2s(1:nact,1:nact,1:nact,1:nact,1:6))
  allocate(den1_noref(1:nact,1:nact,1:2))
  allocate(den2_noref(1:nact,1:nact,1:nact,1:nact,1:6))

  if (cc_rank >= 1) then
     allocate(t1inp((norb1+1):nact,1:norb1,1:2))
     allocate(g1inp(1:norb1,(norb1+1):nact,1:2))
     allocate(dt1inp((norb1+1):nact,1:norb1,1:2))
     allocate(dg1inp(1:norb1,(norb1+1):nact,1:2))
     allocate(t1out((norb1+1):nact,1:norb1,1:2))
     allocate(g1out(1:norb1,(norb1+1):nact,1:2))
     allocate(dt1out((norb1+1):nact,1:norb1,1:2))
     allocate(dg1out(1:norb1,(norb1+1):nact,1:2))
  end if

  if (cc_rank >= 2) then
     allocate(t2inp((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:6))
     allocate(g2inp(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,1:6))
     allocate(dt2inp((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:6))
     allocate(dg2inp(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,1:6))
     allocate(t2out((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:2))
     allocate(g2out(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,1:2))
     allocate(dt2out((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:2))
     allocate(dg2out(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,1:2))
  end if

  if (cc_rank >= 3) then
     allocate(t3inp ((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1,1:20))
     allocate(dt3inp((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1,1:20))
     allocate(g3inp (1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:20))
     allocate(dg3inp(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:20))
     allocate(t3out ((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1,1:2))
     allocate(dt3out((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1,1:2))
     allocate(g3out (1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:2))
     allocate(dg3out(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:2))
!full     allocate(t3inp ((norb1+1):act1_ul,(norb1+1):act1_ul,(norb1+1):act1_ul,&
!full          act1_ll:norb1,act1_ll:norb1,act1_ll:norb1,1:20))
!full     allocate(dt3inp((norb1+1):act1_ul,(norb1+1):act1_ul,(norb1+1):act1_ul,&
!full          act1_ll:norb1,act1_ll:norb1,act1_ll:norb1,1:20))
!full     allocate(g3inp (act1_ll:norb1,act1_ll:norb1,act1_ll:norb1,&
!full          (norb1+1):act1_ul,(norb1+1):act1_ul,(norb1+1):act1_ul,1:20))
!full     allocate(dg3inp(act1_ll:norb1,act1_ll:norb1,act1_ll:norb1,&
!full          (norb1+1):act1_ul,(norb1+1):act1_ul,(norb1+1):act1_ul,1:20))
!full     allocate(t3out ((norb1+1):act1_ul,(norb1+1):act1_ul,(norb1+1):act1_ul,&
!full          act1_ll:norb1,act1_ll:norb1,act1_ll:norb1,1:2))
!full     allocate(dt3out((norb1+1):act1_ul,(norb1+1):act1_ul,(norb1+1):act1_ul,&
!full          act1_ll:norb1,act1_ll:norb1,act1_ll:norb1,1:2))
!full     allocate(g3out (act1_ll:norb1,act1_ll:norb1,act1_ll:norb1,&
!full          (norb1+1):act1_ul,(norb1+1):act1_ul,(norb1+1):act1_ul,1:2))
!full     allocate(dg3out(act1_ll:norb1,act1_ll:norb1,act1_ll:norb1,&
!full          (norb1+1):act1_ul,(norb1+1):act1_ul,(norb1+1):act1_ul,1:2))
  end if

  nXai = 0
  nXij = 0
  nXab = 0
  do irot = 1, nrotaa
     p2 = rotaa_mapb(irot,1)
     p1 = rotaa_mapb(irot,2)
     if (p2 > norb1 .and. p1 <= norb1) then
        nXai = nXai + 1
     else if (p2 >= act1_ll .and. p1 < act1_ll) then
        nXij = nXij + 1
     else if (p2 > act1_ul .and. p1 <= act1_ul) then
        nXab = nXab + 1
     end if
  end do
  if (nXai+nXij+nXab .ne. nrotaa) then
     write(6,"('nXai = ',i5)") nXai
     write(6,"('nXij = ',i5)") nXij
     write(6,"('nXab = ',i5)") nXab
     stop 'tdcc_init: inconsistent nX.'
  end if

  allocate(aiX(nXai))
  allocate(ijX(nXij))
  allocate(abX(nXab))
  nXai = 0
  nXij = 0
  nXab = 0
  do irot = 1, nrotaa
     p2 = rotaa_mapb(irot,1)
     p1 = rotaa_mapb(irot,2)
     if (p2 > norb1 .and. p1 <= norb1) then
        nXai = nXai + 1
        aiX(nXai) = irot
        write(6, "('# aarot_ai:',5i5)") &
             irot,nXai,aiX(nXai),rotaa_mapb(aix(nXai),2),rotaa_mapb(aix(nXai),1)
     else if (p2 >= act1_ll .and. p1 < act1_ll) then
        nXij = nXij + 1
        ijX(nXij) = irot
        write(6, "('# aarot_ij:',5i5)") &
             irot,nXij,ijX(nXij),rotaa_mapb(ijx(nXij),2),rotaa_mapb(ijx(nXij),1)
     else if (p2 > act1_ul .and. p1 <= act1_ul) then
        nXab = nXab + 1
        abX(nXab) = irot
        write(6, "('# aarot_ab:',5i5)") &
             irot,nXab,abX(nXab),rotaa_mapb(abx(nXab),2),rotaa_mapb(abx(nXab),1)
     end if
  end do

  spin_focka = tdcc_spin_fock(1,1)
  spin_t2aa = tdcc_spin_t2inp(1,1,1,1)
  spin_t2ab = tdcc_spin_t2inp(1,2,1,2)
  spin_t3aaa = tdcc_spin_t3inp(1,1,1,1,1,1)
  spin_t3aab = tdcc_spin_t3inp(1,1,2,1,1,2)
  spin_g1a = 1
  spin_g2aa = tdcc_spin_g2inp(1,1,1,1)
  spin_g2ab = tdcc_spin_g2inp(1,2,1,2)
  spin_g3aaa = tdcc_spin_g3inp(1,1,1,1,1,1)
  spin_g3aab = tdcc_spin_g3inp(1,1,2,1,1,2)
  spin_int2aa = tdcc_spin_int2x(1,1,1,1)
  spin_int2ab = tdcc_spin_int2x(1,2,1,2)
  allocate(cc_work1(1:nact**4))
  allocate(cc_work2(1:nact**4))
  allocate(cc_work3(1:nact**4))

end subroutine tdcc_init
!################################################################################
