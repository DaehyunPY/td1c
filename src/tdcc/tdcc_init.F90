!################################################################################
subroutine tdcc_init(mtot)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nelact,nact,ndetx,nsub,norb_sub,max_sub
  use mod_cc, only : cc_type,dot1,optcc,bcc,optbcc
  use mod_cc, only : norb1,norb2,cc_rank,fock,int2x,den1s,den2s,den1_noref,den2_noref
  use mod_cc, only : t1inp,g1inp,t2inp,g2inp,t3inp,g3inp,dt1inp,dg1inp,dt2inp,dg2inp,dt3inp,dg3inp
  use mod_cc, only : t1out,g1out,t2out,g2out,t3out,g3out,dt1out,dg1out,dt2out,dg2out,dt3out,dg3out
  use mod_cc2
  
  implicit none
  integer(c_long), intent(in) :: mtot
  integer,external :: tdcc_spin_int2x,tdcc_spin_fock
  integer,external :: tdcc_spin_t2inp,tdcc_spin_t3inp
  integer,external :: tdcc_spin_g2inp,tdcc_spin_g3inp

  if (ndetx < 2) stop "tdcc_init: ndetx < 2."
  if (nact == 0) stop "tdcc_init: nact = 0."
  if (nsub .ne. 2) stop "tdcc_init: nsub .ne. 2."
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
  end if

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
