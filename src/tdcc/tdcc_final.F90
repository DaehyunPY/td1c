!################################################################################
subroutine tdcc_final()

  use, intrinsic :: iso_c_binding
  use mod_cc, only : cc_rank,map_cc0,map_cc1a,map_cc2aa,map_cc2ab,map_cc3aaa,map_cc3aab
  use mod_cc, only : h1_cc1a,p1_cc1a
  use mod_cc, only : h1_cc2aa,h2_cc2aa,p1_cc2aa,p2_cc2aa
  use mod_cc, only : h1_cc2ab,h2_cc2ab,p1_cc2ab,p2_cc2ab
  use mod_cc, only : h1_cc3aaa,h2_cc3aaa,h3_cc3aaa,p1_cc3aaa,p2_cc3aaa,p3_cc3aaa
  use mod_cc, only : h1_cc3aab,h2_cc3aab,h3_cc3aab,p1_cc3aab,p2_cc3aab,p3_cc3aab
  use mod_cc, only : fock,int2x,den1s,den2s,den1_noref,den2_noref
  use mod_cc, only : t1inp,g1inp,t2inp,g2inp,t3inp,g3inp,dt1inp,dg1inp,dt2inp,dg2inp,dt3inp,dg3inp
  use mod_cc, only : t1out,g1out,t2out,g2out,t3out,g3out,dt1out,dg1out,dt2out,dg2out,dt3out,dg3out
  use mod_cc, only : aiX,ijX,abX
  use mod_cc2
  
  implicit none
  deallocate(map_cc0)
  deallocate(fock)
  deallocate(int2x)
  deallocate(den1s)
  deallocate(den2s)
  deallocate(den1_noref)
  deallocate(den2_noref)

  deallocate(aiX)
  deallocate(ijX)
  deallocate(abX)

  deallocate(cc_work1)
  deallocate(cc_work2)
  deallocate(cc_work3)

  deallocate(h1_ooooaa);deallocate(h2_ooooaa);deallocate(h3_ooooaa);deallocate(h4_ooooaa)
  deallocate(h1_ooooab);deallocate(h2_ooooab);deallocate(h3_ooooab);deallocate(h4_ooooab)
  deallocate(h1_ooovaa);deallocate(h2_ooovaa);deallocate(h3_ooovaa);deallocate(p4_ooovaa)
  deallocate(h1_ooovab);deallocate(h2_ooovab);deallocate(h3_ooovab);deallocate(p4_ooovab)
  deallocate(h1_oovvaa);deallocate(h2_oovvaa);deallocate(p3_oovvaa);deallocate(p4_oovvaa)
  deallocate(h1_oovvab);deallocate(h2_oovvab);deallocate(p3_oovvab);deallocate(p4_oovvab)
  deallocate(h1_ovovaa);deallocate(p2_ovovaa);deallocate(h3_ovovaa);deallocate(p4_ovovaa)
  deallocate(h1_ovovab);deallocate(p2_ovovab);deallocate(h3_ovovab);deallocate(p4_ovovab)
  deallocate(h1_ovvvaa);deallocate(p2_ovvvaa);deallocate(p3_ovvvaa);deallocate(p4_ovvvaa)
  deallocate(h1_ovvvab);deallocate(p2_ovvvab);deallocate(p3_ovvvab);deallocate(p4_ovvvab)
  deallocate(p1_vvvvaa);deallocate(p2_vvvvaa);deallocate(p3_vvvvaa);deallocate(p4_vvvvaa)
  deallocate(p1_vvvvab);deallocate(p2_vvvvab);deallocate(p3_vvvvab);deallocate(p4_vvvvab)

  if (cc_rank >= 1) then
     deallocate(map_cc1a)
     deallocate(h1_cc1a)
     deallocate(p1_cc1a)
     deallocate(t1inp)
     deallocate(g1inp)
     deallocate(dt1inp)
     deallocate(dg1inp)
     deallocate(t1out)
     deallocate(g1out)
     deallocate(dt1out)
     deallocate(dg1out)
  end if
  if (cc_rank >= 2) then
     deallocate(map_cc2aa)
     deallocate(map_cc2ab)
     deallocate(h1_cc2aa)
     deallocate(h2_cc2aa)
     deallocate(p1_cc2aa)
     deallocate(p2_cc2aa)
     deallocate(h1_cc2ab)
     deallocate(h2_cc2ab)
     deallocate(p1_cc2ab)
     deallocate(p2_cc2ab)
     deallocate(t2inp)
     deallocate(g2inp)
     deallocate(dt2inp)
     deallocate(dg2inp)
     deallocate(t2out)
     deallocate(g2out)
     deallocate(dt2out)
     deallocate(dg2out)
  end if
  if (cc_rank >= 3) then
     deallocate(map_cc3aaa)
     deallocate(map_cc3aab)
     deallocate(h1_cc3aaa)
     deallocate(h2_cc3aaa)
     deallocate(h3_cc3aaa)
     deallocate(p1_cc3aaa)
     deallocate(p2_cc3aaa)
     deallocate(p3_cc3aaa)
     deallocate(h1_cc3aab)
     deallocate(h2_cc3aab)
     deallocate(h3_cc3aab)
     deallocate(p1_cc3aab)
     deallocate(p2_cc3aab)
     deallocate(p3_cc3aab)
     deallocate(t3inp)
     deallocate(g3inp)
     deallocate(dt3inp)
     deallocate(dg3inp)
     deallocate(t3out)
     deallocate(g3out)
     deallocate(dt3out)
     deallocate(dg3out)
  end if

end subroutine tdcc_final
!################################################################################
