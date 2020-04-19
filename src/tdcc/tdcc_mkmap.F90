!################################################################################
subroutine tdcc_mkmap()

  use, intrinsic :: iso_c_binding

  use mod_ormas, only : iprint, nact, ncore, nelact, ntot_alph_beta
  use mod_ormas, only : nstr_alph, onv_alph, orb_alph
  use mod_ormas, only : nstr_beta, onv_beta, orb_beta
  use mod_ormas, only : mval_alph, mval_beta
  use mod_ormas, only : llstr_alph_beta, nstr_alph_beta
  use mod_cc, only : cc_rank,ncc0,ncc1a,ncc2aa,ncc2ab,ncc3aaa,ncc3aab,norb1
  use mod_cc, only : map_cc0,map_cc1a,map_cc2aa,map_cc2ab,map_cc3aaa,map_cc3aab
  use mod_cc, only : h1_cc1a,p1_cc1a
  use mod_cc, only : h1_cc2aa,h2_cc2aa,p1_cc2aa,p2_cc2aa
  use mod_cc, only : h1_cc2ab,h2_cc2ab,p1_cc2ab,p2_cc2ab
  use mod_cc, only : h1_cc3aaa,h2_cc3aaa,h3_cc3aaa,p1_cc3aaa,p2_cc3aaa,p3_cc3aaa
  use mod_cc, only : h1_cc3aab,h2_cc3aab,h3_cc3aab,p1_cc3aab,p2_cc3aab,p3_cc3aab

  use mod_ormas, only : act1_ll,act1_ul
  use mod_cc, only : ncc2aa_act1,ncc2ab_act1
  use mod_cc, only : map_cc2aa_act1,map_cc2ab_act1
  use mod_cc, only : h1_cc2aa_act1,h2_cc2aa_act1,p1_cc2aa_act1,p2_cc2aa_act1
  use mod_cc, only : h1_cc2ab_act1,h2_cc2ab_act1,p1_cc2ab_act1,p2_cc2ab_act1

  implicit none
  integer(c_int) :: istr,jstr,ax,bx,idet,icc,a,b,c,d,e,i,j,k,l,m
  integer(c_int) :: ap1,bp1,ah1,bh1

  ! count nonzero tensors
  call tdcc_mkmap_oo
  call tdcc_mkmap_ov
  call tdcc_mkmap_vv
  call tdcc_mkmap_ooooaa
  call tdcc_mkmap_ooooab
  call tdcc_mkmap_ooovaa
  call tdcc_mkmap_ooovab
  call tdcc_mkmap_oovvaa
  call tdcc_mkmap_oovvab
  call tdcc_mkmap_ovovaa
  call tdcc_mkmap_ovovab
  call tdcc_mkmap_ovvvaa
  call tdcc_mkmap_ovvvab
  call tdcc_mkmap_vvvvaa
  call tdcc_mkmap_vvvvab

  ! count nonzero amplitudes
  ncc0 = 0
  ncc1a = 0
  ncc2aa = 0
  ncc2ab = 0
  ncc3aaa = 0
  ncc3aab = 0
  ncc2aa_act1 = 0
  ncc2ab_act1 = 0
  do istr = 1, nstr_beta
     do jstr = llstr_alph_beta(istr), llstr_alph_beta(istr)+nstr_alph_beta(istr)-1
!old        ax = dist_str_alph(1,jstr)-1
!old        bx = dist_str_beta(1,istr)-1
        ax = sum(onv_alph(norb1+1:nact,jstr))
        bx = sum(onv_beta(norb1+1:nact,istr))

        ap1 = sum(onv_alph(norb1+1:act1_ul,jstr))
        bp1 = sum(onv_beta(norb1+1:act1_ul,istr))
        ah1 = norb1-act1_ll+1-sum(onv_alph(act1_ll:norb1,jstr))
        bh1 = norb1-act1_ll+1-sum(onv_beta(act1_ll:norb1,istr))

        idet = ntot_alph_beta(istr)+jstr

        !write(6, "('tdcc: ',3i10,' ')", advance = 'no') idet, jstr, istr
        !call ormas_occvec_print(6, .false., nact,onv_alph(1,jstr))
        !write(6, "(' x ')", advance = 'no')
        !call ormas_occvec_print(6, .false., nact,onv_beta(1,istr))

        if (ax == 0 .and. bx == 0) then
           ncc0 = ncc0 + 1
           !write(6,"('   REF')")
        else if (ax == 1 .and. bx == 0) then
           ncc1a = ncc1a + 1
           !write(6,"('    1a')")
        else if (ax == 2 .and. bx == 0) then
           ncc2aa = ncc2aa + 1
           !write(6,"('   2aa')")
           if (ah1 == 2 .and. ap1 == 2) then
              ncc2aa_act1 = ncc2aa_act1 + 1
           end if
        else if (ax == 1 .and. bx == 1) then
           ncc2ab = ncc2ab + 1
           !write(6,"('   2ab')")
           if (ah1 == 1 .and. ap1 == 1 .and. &
               bh1 == 1 .and. bp1 == 1) then
              ncc2ab_act1 = ncc2ab_act1 + 1
           end if
        else if (ax == 3 .and. bx == 0) then
           ncc3aaa = ncc3aaa + 1
           !write(6,"('  3aaa')")
        else if (ax == 2 .and. bx == 1) then
           ncc3aab = ncc3aab + 1
           !write(6,"('  3aab')")
        else
           !write(6,"('  else')")
        end if

     end do
  end do

  write(6,"('tdcc_mkmap: ncc0    = ',i5)") ncc0
  write(6,"('tdcc_mkmap: ncc1a   = ',i5)") ncc1a
  write(6,"('tdcc_mkmap: ncc2aa  = ',i5)") ncc2aa
  write(6,"('tdcc_mkmap: ncc2ab  = ',i5)") ncc2ab
  write(6,"('tdcc_mkmap: ncc3aaa = ',i5)") ncc3aaa
  write(6,"('tdcc_mkmap: ncc3aab = ',i5)") ncc3aab
  write(6,"('tdcc_mkmap: ncc2aa_act1  = ',i5)") ncc2aa_act1
  write(6,"('tdcc_mkmap: ncc2ab_act1  = ',i5)") ncc2ab_act1
  allocate(map_cc0(1:ncc0,1:3))
  if (cc_rank >= 1) then
     allocate(map_cc1a(1:ncc1a,1:3))
     allocate(h1_cc1a (1:ncc1a))
     allocate(p1_cc1a (1:ncc1a))
  end if
  if (cc_rank >= 2) then
     allocate(map_cc2aa(1:ncc2aa,1:3))
     allocate(map_cc2ab(1:ncc2ab,1:3))
     allocate(h1_cc2aa (1:ncc2aa))
     allocate(h2_cc2aa (1:ncc2aa))
     allocate(p1_cc2aa (1:ncc2aa))
     allocate(p2_cc2aa (1:ncc2aa))
     allocate(h1_cc2ab (1:ncc2ab))
     allocate(h2_cc2ab (1:ncc2ab))
     allocate(p1_cc2ab (1:ncc2ab))
     allocate(p2_cc2ab (1:ncc2ab))
  end if
  if (cc_rank >= 3) then
     allocate(map_cc3aaa(1:ncc3aaa,1:3))
     allocate(map_cc3aab(1:ncc3aab,1:3))
     allocate(h1_cc3aaa (1:ncc3aaa))
     allocate(h2_cc3aaa (1:ncc3aaa))
     allocate(h3_cc3aaa (1:ncc3aaa))
     allocate(p1_cc3aaa (1:ncc3aaa))
     allocate(p2_cc3aaa (1:ncc3aaa))
     allocate(p3_cc3aaa (1:ncc3aaa))
     allocate(h1_cc3aab (1:ncc3aab))
     allocate(h2_cc3aab (1:ncc3aab))
     allocate(h3_cc3aab (1:ncc3aab))
     allocate(p1_cc3aab (1:ncc3aab))
     allocate(p2_cc3aab (1:ncc3aab))
     allocate(p3_cc3aab (1:ncc3aab))
     allocate(map_cc2aa_act1(1:ncc2aa_act1,1:3))
     allocate(map_cc2ab_act1(1:ncc2ab_act1,1:3))
     allocate(h1_cc2aa_act1 (1:ncc2aa_act1))
     allocate(h2_cc2aa_act1 (1:ncc2aa_act1))
     allocate(p1_cc2aa_act1 (1:ncc2aa_act1))
     allocate(p2_cc2aa_act1 (1:ncc2aa_act1))
     allocate(h1_cc2ab_act1 (1:ncc2ab_act1))
     allocate(h2_cc2ab_act1 (1:ncc2ab_act1))
     allocate(p1_cc2ab_act1 (1:ncc2ab_act1))
     allocate(p2_cc2ab_act1 (1:ncc2ab_act1))
  end if

  ncc0 = 0
  ncc1a = 0
  ncc2aa = 0
  ncc2ab = 0
  ncc3aaa = 0
  ncc3aab = 0
  ncc2aa_act1 = 0
  ncc2ab_act1 = 0

  do istr = 1, nstr_beta
     do jstr = llstr_alph_beta(istr), llstr_alph_beta(istr)+nstr_alph_beta(istr)-1
!old        ax = dist_str_alph(1,jstr)-1
!old        bx = dist_str_beta(1,istr)-1
        ax = sum(onv_alph(norb1+1:nact,jstr))
        bx = sum(onv_beta(norb1+1:nact,istr))
        ap1 = sum(onv_alph(norb1+1:act1_ul,jstr))
        bp1 = sum(onv_beta(norb1+1:act1_ul,istr))
        ah1 = norb1-act1_ll+1-sum(onv_alph(act1_ll:norb1,jstr))
        bh1 = norb1-act1_ll+1-sum(onv_beta(act1_ll:norb1,istr))

        idet = ntot_alph_beta(istr)+jstr
        if (ax == 0 .and. bx == 0) then
           ncc0 = ncc0 + 1
           map_cc0(ncc0,1) = jstr
           map_cc0(ncc0,2) = istr
           map_cc0(ncc0,3) = idet
        else if (ax == 1 .and. bx == 0) then
           ncc1a = ncc1a + 1
           map_cc1a(ncc1a,1) = jstr
           map_cc1a(ncc1a,2) = istr
           map_cc1a(ncc1a,3) = idet
           call tdcc_mkmap_hpcc1a(onv_alph(1,jstr), &
                h1_cc1a(ncc1a), &
                p1_cc1a(ncc1a))
        else if (ax == 2 .and. bx == 0) then
           ncc2aa = ncc2aa + 1
           map_cc2aa(ncc2aa,1) = jstr
           map_cc2aa(ncc2aa,2) = istr
           map_cc2aa(ncc2aa,3) = idet
           call tdcc_mkmap_hpcc2aa(onv_alph(1,jstr), &
                h1_cc2aa(ncc2aa),h2_cc2aa(ncc2aa), &
                p1_cc2aa(ncc2aa),p2_cc2aa(ncc2aa))
           if (cc_rank >= 3 .and. &
               ah1 == 2 .and. ap1 == 2) then
              ncc2aa_act1 = ncc2aa_act1 + 1
              map_cc2aa_act1(ncc2aa_act1,1) = jstr
              map_cc2aa_act1(ncc2aa_act1,2) = istr
              map_cc2aa_act1(ncc2aa_act1,3) = idet
              call tdcc_mkmap_hpcc2aa_act1(onv_alph(1,jstr), &
                   h1_cc2aa_act1(ncc2aa_act1),h2_cc2aa_act1(ncc2aa_act1), &
                   p1_cc2aa_act1(ncc2aa_act1),p2_cc2aa_act1(ncc2aa_act1))
           end if
        else if (ax == 1 .and. bx == 1) then
           ncc2ab = ncc2ab + 1
           map_cc2ab(ncc2ab,1) = jstr
           map_cc2ab(ncc2ab,2) = istr
           map_cc2ab(ncc2ab,3) = idet
           call tdcc_mkmap_hpcc2ab(onv_alph(1,jstr),onv_beta(1,istr), &
                h1_cc2ab(ncc2ab),h2_cc2ab(ncc2ab), &
                p1_cc2ab(ncc2ab),p2_cc2ab(ncc2ab))
           if (cc_rank >= 3 .and. &
               ah1 == 1 .and. ap1 == 1 .and. &
               bh1 == 1 .and. bp1 == 1) then
              ncc2ab_act1 = ncc2ab_act1 + 1
              map_cc2ab_act1(ncc2ab_act1,1) = jstr
              map_cc2ab_act1(ncc2ab_act1,2) = istr
              map_cc2ab_act1(ncc2ab_act1,3) = idet
              call tdcc_mkmap_hpcc2ab_act1(onv_alph(1,jstr),onv_beta(1,istr), &
                   h1_cc2ab_act1(ncc2ab_act1),h2_cc2ab_act1(ncc2ab_act1), &
                   p1_cc2ab_act1(ncc2ab_act1),p2_cc2ab_act1(ncc2ab_act1))
           end if
        else if (ax == 3 .and. bx == 0) then
           ncc3aaa = ncc3aaa + 1
           map_cc3aaa(ncc3aaa,1) = jstr
           map_cc3aaa(ncc3aaa,2) = istr
           map_cc3aaa(ncc3aaa,3) = idet
           call tdcc_mkmap_hpcc3aaa(onv_alph(1,jstr), &
                h1_cc3aaa(ncc3aaa),h2_cc3aaa(ncc3aaa),h3_cc3aaa(ncc3aaa), &
                p1_cc3aaa(ncc3aaa),p2_cc3aaa(ncc3aaa),p3_cc3aaa(ncc3aaa))
        else if (ax == 2 .and. bx == 1) then
           ncc3aab = ncc3aab + 1
           map_cc3aab(ncc3aab,1) = jstr
           map_cc3aab(ncc3aab,2) = istr
           map_cc3aab(ncc3aab,3) = idet
           call tdcc_mkmap_hpcc3aab(onv_alph(1,jstr),onv_beta(1,istr), &
                h1_cc3aab(ncc3aab),h2_cc3aab(ncc3aab),h3_cc3aab(ncc3aab), &
                p1_cc3aab(ncc3aab),p2_cc3aab(ncc3aab),p3_cc3aab(ncc3aab))
        end if
     end do
  end do

  if (iprint > 1) then
     do icc = 1, ncc0
        write(6,"('tdcc_mkmap:  REF ',4i10)",advance='no') icc,map_cc0(icc,1),map_cc0(icc,2),map_cc0(icc,3)
        write(6, "('   ')", advance = 'no')
        call ormas_occvec_print(6, .false., nact,onv_alph(1,map_cc0(icc,1)))
        write(6, "(' x ')", advance = 'no')
        call ormas_occvec_print(6, .true., nact,onv_beta(1,map_cc0(icc,2)))
     end do
     do icc = 1, ncc1a
        write(6,"('tdcc_mkmap:   1a ',4i10)",advance='no') icc,map_cc1a(icc,1),map_cc1a(icc,2),map_cc1a(icc,3)
        write(6, "('   ')", advance = 'no')
        call ormas_occvec_print(6, .false., nact,onv_alph(1,map_cc1a(icc,1)))
        write(6, "(' x ')", advance = 'no')
        call ormas_occvec_print(6, .false., nact,onv_beta(1,map_cc1a(icc,2)))
        write(6,"(2i5)") h1_cc1a(icc), &
                         p1_cc1a(icc)
     end do
     do icc = 1, ncc2aa
        write(6,"('tdcc_mkmap:  2aa ',4i10)",advance='no') icc,map_cc2aa(icc,1),map_cc2aa(icc,2),map_cc2aa(icc,3)
        write(6, "('   ')", advance = 'no')
        call ormas_occvec_print(6, .false., nact,onv_alph(1,map_cc2aa(icc,1)))
        write(6, "(' x ')", advance = 'no')
        call ormas_occvec_print(6, .false., nact,onv_beta(1,map_cc2aa(icc,2)))
        write(6,"(4i5)") h1_cc2aa(icc),h2_cc2aa(icc), &
                         p1_cc2aa(icc),p2_cc2aa(icc)
     end do
     do icc = 1, ncc2ab
        write(6,"('tdcc_mkmap:  2ab ',4i10)",advance='no') icc,map_cc2ab(icc,1),map_cc2ab(icc,2),map_cc2ab(icc,3)
        write(6, "('   ')", advance = 'no')
        call ormas_occvec_print(6, .false., nact,onv_alph(1,map_cc2ab(icc,1)))
        write(6, "(' x ')", advance = 'no')
        call ormas_occvec_print(6, .false., nact,onv_beta(1,map_cc2ab(icc,2)))
        write(6,"(4i5)") h1_cc2ab(icc),h2_cc2ab(icc), &
                         p1_cc2ab(icc),p2_cc2ab(icc)
     end do
     if (cc_rank >= 3) then
        do icc = 1, ncc2aa_act1
           write(6,"('tdcc_mkmap:  2aa_act1 ',4i10)",advance='no') &
                icc,map_cc2aa_act1(icc,1),map_cc2aa_act1(icc,2),map_cc2aa_act1(icc,3)
           write(6, "('   ')", advance = 'no')
           call ormas_occvec_print(6, .false., nact,onv_alph(1,map_cc2aa_act1(icc,1)))
           write(6, "(' x ')", advance = 'no')
           call ormas_occvec_print(6, .false., nact,onv_beta(1,map_cc2aa_act1(icc,2)))
           write(6,"(4i5)") h1_cc2aa_act1(icc),h2_cc2aa_act1(icc), &
                            p1_cc2aa_act1(icc),p2_cc2aa_act1(icc)
        end do
        do icc = 1, ncc2ab_act1
           write(6,"('tdcc_mkmap:  2ab_act1 ',4i10)",advance='no') &
                icc,map_cc2ab_act1(icc,1),map_cc2ab_act1(icc,2),map_cc2ab_act1(icc,3)
           write(6, "('   ')", advance = 'no')
           call ormas_occvec_print(6, .false., nact,onv_alph(1,map_cc2ab_act1(icc,1)))
           write(6, "(' x ')", advance = 'no')
           call ormas_occvec_print(6, .false., nact,onv_beta(1,map_cc2ab_act1(icc,2)))
           write(6,"(4i5)") h1_cc2ab_act1(icc),h2_cc2ab_act1(icc), &
                            p1_cc2ab_act1(icc),p2_cc2ab_act1(icc)
        end do
     end if
     do icc = 1, ncc3aaa
        write(6,"('tdcc_mkmap: 3aaa ',4i10)",advance='no') icc,map_cc3aaa(icc,1),map_cc3aaa(icc,2),map_cc3aaa(icc,3)
        write(6, "('   ')", advance = 'no')
        call ormas_occvec_print(6, .false., nact,onv_alph(1,map_cc3aaa(icc,1)))
        write(6, "(' x ')", advance = 'no')
        call ormas_occvec_print(6, .false., nact,onv_beta(1,map_cc3aaa(icc,2)))
        write(6,"(6i5)") h1_cc3aaa(icc),h2_cc3aaa(icc),h3_cc3aaa(icc), &
                         p1_cc3aaa(icc),p2_cc3aaa(icc),p3_cc3aaa(icc)
     end do
     do icc = 1, ncc3aab
        write(6,"('tdcc_mkmap: 3aab ',4i10)",advance='no') icc,map_cc3aab(icc,1),map_cc3aab(icc,2),map_cc3aab(icc,3)
        write(6, "('   ')", advance = 'no')
        call ormas_occvec_print(6, .false., nact,onv_alph(1,map_cc3aab(icc,1)))
        write(6, "(' x ')", advance = 'no')
        call ormas_occvec_print(6, .false., nact,onv_beta(1,map_cc3aab(icc,2)))
        write(6,"(6i5)") h1_cc3aab(icc),h2_cc3aab(icc),h3_cc3aab(icc), &
                         p1_cc3aab(icc),p2_cc3aab(icc),p3_cc3aab(icc)
     end do
  end if

end subroutine tdcc_mkmap
!################################################################################
subroutine tdcc_mkmap_hpcc1a(onva,h1,p1)
  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nact
  use mod_cc, only : norb1
  use mod_cc2

  implicit none
  integer(c_int), intent(in) :: onva(1:nact)
  integer(c_int), intent(out) :: h1,p1
  integer(c_int) :: iact

  do iact = 1, norb1
     if (onva(iact) == 0) then
        h1 = iact
        exit
     end if
  end do

  do iact = norb1+1, nact
     if (onva(iact) == 1) then
        p1 = iact
        exit
     end if
  end do

end subroutine tdcc_mkmap_hpcc1a
!################################################################################
subroutine tdcc_mkmap_hpcc2aa(onva,h1,h2,p1,p2)
  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nact
  use mod_cc, only : norb1
  use mod_cc2

  implicit none
  integer(c_int), intent(in) :: onva(1:nact)
  integer(c_int), intent(out) :: h1,h2,p1,p2
  integer(c_int) :: iact

  do iact = 1, norb1
     if (onva(iact) == 0) then
        h1 = iact
        exit
     end if
  end do
  do iact = h1+1, norb1
     if (onva(iact) == 0) then
        h2 = iact
        exit
     end if
  end do

  do iact = norb1+1, nact
     if (onva(iact) == 1) then
        p1 = iact
        exit
     end if
  end do
  do iact = p1+1, nact
     if (onva(iact) == 1) then
        p2 = iact
        exit
     end if
  end do

end subroutine tdcc_mkmap_hpcc2aa
!################################################################################
subroutine tdcc_mkmap_hpcc2ab(onva,onvb,h1,h2,p1,p2)
  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nact
  use mod_cc, only : norb1
  use mod_cc2

  implicit none
  integer(c_int), intent(in) :: onva(1:nact),onvb(1:nact)
  integer(c_int), intent(out) :: h1,h2,p1,p2
  integer(c_int) :: iact

  do iact = 1, norb1
     if (onva(iact) == 0) then
        h1 = iact
        exit
     end if
  end do
  do iact = 1, norb1
     if (onvb(iact) == 0) then
        h2 = iact
        exit
     end if
  end do

  do iact = norb1+1, nact
     if (onva(iact) == 1) then
        p1 = iact
        exit
     end if
  end do
  do iact = norb1+1, nact
     if (onvb(iact) == 1) then
        p2 = iact
        exit
     end if
  end do

end subroutine tdcc_mkmap_hpcc2ab
!################################################################################
subroutine tdcc_mkmap_hpcc2aa_act1(onva,h1,h2,p1,p2)
  use, intrinsic :: iso_c_binding
  use mod_ormas, only : act1_ll,act1_ul
  use mod_ormas, only : nact
  use mod_cc, only : norb1
  use mod_cc2

  implicit none
  integer(c_int), intent(in) :: onva(1:nact)
  integer(c_int), intent(out) :: h1,h2,p1,p2
  integer(c_int) :: iact

  do iact = act1_ll, norb1
     if (onva(iact) == 0) then
        h1 = iact
        exit
     end if
  end do
  do iact = h1+1, norb1
     if (onva(iact) == 0) then
        h2 = iact
        exit
     end if
  end do

  do iact = norb1+1, act1_ul
     if (onva(iact) == 1) then
        p1 = iact
        exit
     end if
  end do
  do iact = p1+1, act1_ul
     if (onva(iact) == 1) then
        p2 = iact
        exit
     end if
  end do

end subroutine tdcc_mkmap_hpcc2aa_act1
!################################################################################
subroutine tdcc_mkmap_hpcc2ab_act1(onva,onvb,h1,h2,p1,p2)
  use, intrinsic :: iso_c_binding
  use mod_ormas, only : act1_ll,act1_ul
  use mod_ormas, only : nact
  use mod_cc, only : norb1
  use mod_cc2

  implicit none
  integer(c_int), intent(in) :: onva(1:nact),onvb(1:nact)
  integer(c_int), intent(out) :: h1,h2,p1,p2
  integer(c_int) :: iact

  do iact = act1_ll, norb1
     if (onva(iact) == 0) then
        h1 = iact
        exit
     end if
  end do
  do iact = act1_ll, norb1
     if (onvb(iact) == 0) then
        h2 = iact
        exit
     end if
  end do

  do iact = norb1+1, act1_ul
     if (onva(iact) == 1) then
        p1 = iact
        exit
     end if
  end do
  do iact = norb1+1, act1_ul
     if (onvb(iact) == 1) then
        p2 = iact
        exit
     end if
  end do

end subroutine tdcc_mkmap_hpcc2ab_act1
!################################################################################
subroutine tdcc_mkmap_hpcc3aaa(onva,h1,h2,h3,p1,p2,p3)
  use, intrinsic :: iso_c_binding
  use mod_ormas, only : act1_ll,act1_ul
  use mod_cc, only : norb1
  use mod_cc2

  implicit none
  integer(c_int), intent(in) :: onva(1:*)
  integer(c_int), intent(out) :: h1,h2,h3,p1,p2,p3
  integer(c_int) :: iact

!dplus
!  do iact = 1, norb1
  do iact = act1_ll, norb1
     if (onva(iact) == 0) then
        h1 = iact
        exit
     end if
  end do
  do iact = h1+1, norb1
     if (onva(iact) == 0) then
        h2 = iact
        exit
     end if
  end do
  do iact = h2+1, norb1
     if (onva(iact) == 0) then
        h3 = iact
        exit
     end if
  end do

!dplus
!  do iact = norb1+1, nact
  do iact = norb1+1, act1_ul
     if (onva(iact) == 1) then
        p1 = iact
        exit
     end if
  end do
!dplus
!  do iact = p1+1, nact
  do iact = p1+1, act1_ul
     if (onva(iact) == 1) then
        p2 = iact
        exit
     end if
  end do
!dplus
!  do iact = p2+1, nact
  do iact = p2+1, act1_ul
     if (onva(iact) == 1) then
        p3 = iact
        exit
     end if
  end do

end subroutine tdcc_mkmap_hpcc3aaa
!################################################################################
subroutine tdcc_mkmap_hpcc3aab(onva,onvb,h1,h2,h3,p1,p2,p3)
  use, intrinsic :: iso_c_binding
  use mod_ormas, only : act1_ll,act1_ul
  use mod_cc, only : norb1
  use mod_cc2

  implicit none
  integer(c_int), intent(in) :: onva(1:*),onvb(1:*)
  integer(c_int), intent(out) :: h1,h2,h3,p1,p2,p3
  integer(c_int) :: iact

!dplus
!  do iact = 1, norb1
  do iact = act1_ll, norb1
     if (onva(iact) == 0) then
        h1 = iact
        exit
     end if
  end do
  do iact = h1+1, norb1
     if (onva(iact) == 0) then
        h2 = iact
        exit
     end if
  end do
!dplus
!  do iact = 1, norb1
  do iact = act1_ll, norb1
     if (onvb(iact) == 0) then
        h3 = iact
        exit
     end if
  end do

!dplus
!  do iact = norb1+1, nact
  do iact = norb1+1, act1_ul
     if (onva(iact) == 1) then
        p1 = iact
        exit
     end if
  end do
!dplus
!  do iact = p1+1, nact
  do iact = p1+1, act1_ul
     if (onva(iact) == 1) then
        p2 = iact
        exit
     end if
  end do
!dplus
!  do iact = norb1+1, nact
  do iact = norb1+1, act1_ul
     if (onvb(iact) == 1) then
        p3 = iact
        exit
     end if
  end do

end subroutine tdcc_mkmap_hpcc3aab
!################################################################################
subroutine tdcc_mkmap_oo()

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : mval,nact,ncore
  use mod_cc, only : norb1
  use mod_cc2

  implicit none
  integer(c_int) :: i,j

  noo = 0
  do i = 1,norb1
  do j = 1,norb1
     if (mval(ncore+i)==mval(ncore+j)) then
        noo = noo + 1
     end if
  end do
  end do
  write(6,"('tdcc_mkmap: noo     = ',2i5)") noo,norb1**2
  allocate(h1_oo(1:noo))
  allocate(h2_oo(1:noo))
  noo = 0
  do i = 1,norb1
  do j = 1,norb1
     if (mval(ncore+i)==mval(ncore+j)) then
        noo = noo + 1
        h1_oo(noo) = i
        h2_oo(noo) = j
     end if
  end do
  end do

end subroutine tdcc_mkmap_oo
!################################################################################
subroutine tdcc_mkmap_ov()

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : mval,nact,ncore
  use mod_cc, only : norb1,norb2
  use mod_cc2

  implicit none
  integer(c_int) :: i,a

  nov = 0
  do i = 1,norb1
  do a = norb1+1,nact
     if (mval(ncore+i)==mval(ncore+a)) then
        nov = nov + 1
     end if
  end do
  end do
  write(6,"('tdcc_mkmap: nov     = ',2i5)") nov,norb1*norb2
  allocate(h1_ov(1:nov))
  allocate(p2_ov(1:nov))
  nov = 0
  do i = 1,norb1
  do a = norb1+1,nact
     if (mval(ncore+i)==mval(ncore+a)) then
        nov = nov + 1
        h1_ov(nov) = i
        p2_ov(nov) = a
     end if
  end do
  end do

end subroutine tdcc_mkmap_ov
!################################################################################
subroutine tdcc_mkmap_vv()

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : mval,nact,ncore
  use mod_cc, only : norb1,norb2
  use mod_cc2

  implicit none
  integer(c_int) :: a,b

  nvv = 0
  do a = norb1+1,nact
  do b = norb1+1,nact
     if (mval(ncore+a)==mval(ncore+b)) then
        nvv = nvv + 1
     end if
  end do
  end do
  write(6,"('tdcc_mkmap: nvv     = ',2i5)") nvv,norb2**2
  allocate(p1_vv(1:nvv))
  allocate(p2_vv(1:nvv))
  nvv = 0
  do a = norb1+1,nact
  do b = norb1+1,nact
     if (mval(ncore+a)==mval(ncore+b)) then
        nvv = nvv + 1
        p1_vv(nvv) = a
        p2_vv(nvv) = b
     end if
  end do
  end do

end subroutine tdcc_mkmap_vv
!################################################################################
subroutine tdcc_mkmap_ooooaa()

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : mval,nact,ncore
  use mod_cc, only : norb1
  use mod_cc2

  implicit none
  integer(c_int) :: i,j,k,l

  ! count nonzero tensors
  nooooaa = 0
  do i = 1,norb1
  do j = 1,i-1
     do k = 1,norb1
     do l = 1,k-1
        if (mval(ncore+i)+mval(ncore+j)==mval(ncore+k)+mval(ncore+l)) then
           nooooaa = nooooaa + 1
        end if
     end do
     end do
  end do
  end do
  write(6,"('tdcc_mkmap: nooooaa = ',2i5)") nooooaa,norb1**4

  allocate(h1_ooooaa(1:nooooaa))
  allocate(h2_ooooaa(1:nooooaa))
  allocate(h3_ooooaa(1:nooooaa))
  allocate(h4_ooooaa(1:nooooaa))

  nooooaa = 0
  do i = 1,norb1
  do j = 1,i-1
     do k = 1,norb1
     do l = 1,k-1
        if (mval(ncore+i)+mval(ncore+j)==mval(ncore+k)+mval(ncore+l)) then
           nooooaa = nooooaa + 1
           h1_ooooaa(nooooaa) = i
           h2_ooooaa(nooooaa) = j
           h3_ooooaa(nooooaa) = k
           h4_ooooaa(nooooaa) = l
        end if
     end do
     end do
  end do
  end do

end subroutine tdcc_mkmap_ooooaa
!################################################################################
subroutine tdcc_mkmap_ooooab()

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : mval,nact,ncore
  use mod_cc, only : norb1
  use mod_cc2

  implicit none
  integer(c_int) :: i,j,k,l

  ! count nonzero tensors
  nooooab = 0
  do i = 1,norb1
  do j = 1,norb1
     do k = 1,norb1
     do l = 1,norb1
        if (mval(ncore+i)+mval(ncore+j)==mval(ncore+k)+mval(ncore+l)) then
           nooooab = nooooab + 1
        end if
     end do
     end do
  end do
  end do
  write(6,"('tdcc_mkmap: nooooab = ',2i5)") nooooab,norb1**4

  allocate(h1_ooooab(1:nooooab))
  allocate(h2_ooooab(1:nooooab))
  allocate(h3_ooooab(1:nooooab))
  allocate(h4_ooooab(1:nooooab))

  nooooab = 0
  do i = 1,norb1
  do j = 1,norb1
     do k = 1,norb1
     do l = 1,norb1
        if (mval(ncore+i)+mval(ncore+j)==mval(ncore+k)+mval(ncore+l)) then
           nooooab = nooooab + 1
           h1_ooooab(nooooab) = i
           h2_ooooab(nooooab) = j
           h3_ooooab(nooooab) = k
           h4_ooooab(nooooab) = l
        end if
     end do
     end do
  end do
  end do

end subroutine tdcc_mkmap_ooooab
!################################################################################
subroutine tdcc_mkmap_ooovaa()

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : mval,nact,ncore
  use mod_cc, only : norb1,norb2
  use mod_cc2

  implicit none
  integer(c_int) :: i,j,k,a

  ! count nonzero tensors
  nooovaa = 0
  do i = 1,norb1
  do j = 1,i-1
     do k = 1,norb1
     do a = norb1+1,nact
        if (mval(ncore+i)+mval(ncore+j)==mval(ncore+k)+mval(ncore+a)) then
           nooovaa = nooovaa + 1
        end if
     end do
     end do
  end do
  end do
  write(6,"('tdcc_mkmap: nooovaa = ',2i5)") nooovaa,norb1**3*norb2

  allocate(h1_ooovaa(1:nooovaa))
  allocate(h2_ooovaa(1:nooovaa))
  allocate(h3_ooovaa(1:nooovaa))
  allocate(p4_ooovaa(1:nooovaa))

  nooovaa = 0
  do i = 1,norb1
  do j = 1,i-1
     do k = 1,norb1
     do a = norb1+1,nact
        if (mval(ncore+i)+mval(ncore+j)==mval(ncore+k)+mval(ncore+a)) then
           nooovaa = nooovaa + 1
           h1_ooovaa(nooovaa) = i
           h2_ooovaa(nooovaa) = j
           h3_ooovaa(nooovaa) = k
           p4_ooovaa(nooovaa) = a
        end if
     end do
     end do
  end do
  end do

end subroutine tdcc_mkmap_ooovaa
!################################################################################
subroutine tdcc_mkmap_ooovab()

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : mval,nact,ncore
  use mod_cc, only : norb1,norb2
  use mod_cc2

  implicit none
  integer(c_int) :: i,j,k,a

  ! count nonzero tensors
  nooovab = 0
  do i = 1,norb1
  do j = 1,norb1
     do k = 1,norb1
     do a = norb1+1,nact
        if (mval(ncore+i)+mval(ncore+j)==mval(ncore+k)+mval(ncore+a)) then
           nooovab = nooovab + 1
        end if
     end do
     end do
  end do
  end do
  write(6,"('tdcc_mkmap: nooovab = ',2i5)") nooovab,norb1**3*norb2

  allocate(h1_ooovab(1:nooovab))
  allocate(h2_ooovab(1:nooovab))
  allocate(h3_ooovab(1:nooovab))
  allocate(p4_ooovab(1:nooovab))

  nooovab = 0
  do i = 1,norb1
  do j = 1,norb1
     do k = 1,norb1
     do a = norb1+1,nact
        if (mval(ncore+i)+mval(ncore+j)==mval(ncore+k)+mval(ncore+a)) then
           nooovab = nooovab + 1
           h1_ooovab(nooovab) = i
           h2_ooovab(nooovab) = j
           h3_ooovab(nooovab) = k
           p4_ooovab(nooovab) = a
        end if
     end do
     end do
  end do
  end do

end subroutine tdcc_mkmap_ooovab
!################################################################################
subroutine tdcc_mkmap_oovvaa()

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : mval,nact,ncore
  use mod_cc, only : norb1,norb2
  use mod_cc2

  implicit none
  integer(c_int) :: i,j,a,b

  ! count nonzero tensors
  noovvaa = 0
  do i = 1,norb1
  do j = 1,i-1
     do a = norb1+1,nact
     do b = norb1+1,a-1
        if (mval(ncore+i)+mval(ncore+j)==mval(ncore+a)+mval(ncore+b)) then
           noovvaa = noovvaa + 1
        end if
     end do
     end do
  end do
  end do
  write(6,"('tdcc_mkmap: noovvaa = ',2i5)") noovvaa,norb1**2*norb2**2

  allocate(h1_oovvaa(1:noovvaa))
  allocate(h2_oovvaa(1:noovvaa))
  allocate(p3_oovvaa(1:noovvaa))
  allocate(p4_oovvaa(1:noovvaa))

  noovvaa = 0
  do i = 1,norb1
  do j = 1,i-1
     do a = norb1+1,nact
     do b = norb1+1,a-1
        if (mval(ncore+i)+mval(ncore+j)==mval(ncore+a)+mval(ncore+b)) then
           noovvaa = noovvaa + 1
           h1_oovvaa(noovvaa) = i
           h2_oovvaa(noovvaa) = j
           p3_oovvaa(noovvaa) = a
           p4_oovvaa(noovvaa) = b
        end if
     end do
     end do
  end do
  end do

end subroutine tdcc_mkmap_oovvaa
!################################################################################
subroutine tdcc_mkmap_oovvab()

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : mval,nact,ncore
  use mod_cc, only : norb1,norb2
  use mod_cc2

  implicit none
  integer(c_int) :: i,j,a,b

  ! count nonzero tensors
  noovvab = 0
  do i = 1,norb1
  do j = 1,norb1
     do a = norb1+1,nact
     do b = norb1+1,nact
        if (mval(ncore+i)+mval(ncore+j)==mval(ncore+a)+mval(ncore+b)) then
           noovvab = noovvab + 1
        end if
     end do
     end do
  end do
  end do
  write(6,"('tdcc_mkmap: noovvab = ',2i5)") noovvab,norb1**2*norb2**2

  allocate(h1_oovvab(1:noovvab))
  allocate(h2_oovvab(1:noovvab))
  allocate(p3_oovvab(1:noovvab))
  allocate(p4_oovvab(1:noovvab))

  noovvab = 0
  do i = 1,norb1
  do j = 1,norb1
     do a = norb1+1,nact
     do b = norb1+1,nact
        if (mval(ncore+i)+mval(ncore+j)==mval(ncore+a)+mval(ncore+b)) then
           noovvab = noovvab + 1
           h1_oovvab(noovvab) = i
           h2_oovvab(noovvab) = j
           p3_oovvab(noovvab) = a
           p4_oovvab(noovvab) = b
        end if
     end do
     end do
  end do
  end do

end subroutine tdcc_mkmap_oovvab
!################################################################################
subroutine tdcc_mkmap_ovovaa()

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : mval,nact,ncore
  use mod_cc, only : norb1,norb2
  use mod_cc2

  implicit none
  integer(c_int) :: i,j,a,b

  ! count nonzero tensors
  novovaa = 0
  do i = 1,norb1
  do a = norb1+1,nact
     do j = 1,norb1
     do b = norb1+1,nact
        if (mval(ncore+i)+mval(ncore+a)==mval(ncore+j)+mval(ncore+b)) then
           novovaa = novovaa + 1
        end if
     end do
     end do
  end do
  end do
  write(6,"('tdcc_mkmap: novovaa = ',2i5)") novovaa,norb1**2*norb2**2

  allocate(h1_ovovaa(1:novovaa))
  allocate(p2_ovovaa(1:novovaa))
  allocate(h3_ovovaa(1:novovaa))
  allocate(p4_ovovaa(1:novovaa))

  novovaa = 0
  do i = 1,norb1
  do a = norb1+1,nact
     do j = 1,norb1
     do b = norb1+1,nact
        if (mval(ncore+i)+mval(ncore+a)==mval(ncore+j)+mval(ncore+b)) then
           novovaa = novovaa + 1
           h1_ovovaa(novovaa) = i
           p2_ovovaa(novovaa) = a
           h3_ovovaa(novovaa) = j
           p4_ovovaa(novovaa) = b
        end if
     end do
     end do
  end do
  end do

end subroutine tdcc_mkmap_ovovaa
!################################################################################
subroutine tdcc_mkmap_ovovab()

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : mval,nact,ncore
  use mod_cc, only : norb1,norb2
  use mod_cc2

  implicit none
  integer(c_int) :: i,j,a,b

  ! count nonzero tensors
  novovab = 0
  do i = 1,norb1
  do a = norb1+1,nact
     do j = 1,norb1
     do b = norb1+1,nact
        if (mval(ncore+i)+mval(ncore+a)==mval(ncore+j)+mval(ncore+b)) then
           novovab = novovab + 1
        end if
     end do
     end do
  end do
  end do
  write(6,"('tdcc_mkmap: novovab = ',2i5)") novovab,norb1**2*norb2**2

  allocate(h1_ovovab(1:novovab))
  allocate(p2_ovovab(1:novovab))
  allocate(h3_ovovab(1:novovab))
  allocate(p4_ovovab(1:novovab))

  novovab = 0
  do i = 1,norb1
  do a = norb1+1,nact
     do j = 1,norb1
     do b = norb1+1,nact
        if (mval(ncore+i)+mval(ncore+a)==mval(ncore+j)+mval(ncore+b)) then
           novovab = novovab + 1
           h1_ovovab(novovab) = i
           p2_ovovab(novovab) = a
           h3_ovovab(novovab) = j
           p4_ovovab(novovab) = b
        end if
     end do
     end do
  end do
  end do

end subroutine tdcc_mkmap_ovovab
!################################################################################
subroutine tdcc_mkmap_ovvvaa()

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : mval,nact,ncore
  use mod_cc, only : norb1,norb2
  use mod_cc2

  implicit none
  integer(c_int) :: i,b,c,d

  ! count nonzero tensors
  novvvaa = 0
  do i = 1,norb1
  do b = norb1+1,nact
     do c = norb1+1,nact
     do d = norb1+1,c-1
        if (mval(ncore+i)+mval(ncore+b)==mval(ncore+c)+mval(ncore+d)) then
           novvvaa = novvvaa + 1
        end if
     end do
     end do
  end do
  end do
  write(6,"('tdcc_mkmap: novvvaa = ',2i5)") novvvaa,norb1*norb2**3

  allocate(h1_ovvvaa(1:novvvaa))
  allocate(p2_ovvvaa(1:novvvaa))
  allocate(p3_ovvvaa(1:novvvaa))
  allocate(p4_ovvvaa(1:novvvaa))

  novvvaa = 0
  do i = 1,norb1
  do b = norb1+1,nact
     do c = norb1+1,nact
     do d = norb1+1,c-1
        if (mval(ncore+i)+mval(ncore+b)==mval(ncore+c)+mval(ncore+d)) then
           novvvaa = novvvaa + 1
           h1_ovvvaa(novvvaa) = i
           p2_ovvvaa(novvvaa) = b
           p3_ovvvaa(novvvaa) = c
           p4_ovvvaa(novvvaa) = d
        end if
     end do
     end do
  end do
  end do

end subroutine tdcc_mkmap_ovvvaa
!################################################################################
subroutine tdcc_mkmap_ovvvab()

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : mval,nact,ncore
  use mod_cc, only : norb1,norb2
  use mod_cc2

  implicit none
  integer(c_int) :: i,b,c,d

  ! count nonzero tensors
  novvvab = 0
  do i = 1,norb1
  do b = norb1+1,nact
     do c = norb1+1,nact
     do d = norb1+1,nact
        if (mval(ncore+i)+mval(ncore+b)==mval(ncore+c)+mval(ncore+d)) then
           novvvab = novvvab + 1
        end if
     end do
     end do
  end do
  end do
  write(6,"('tdcc_mkmap: novvvab = ',2i5)") novvvab,norb1*norb2**3

  allocate(h1_ovvvab(1:novvvab))
  allocate(p2_ovvvab(1:novvvab))
  allocate(p3_ovvvab(1:novvvab))
  allocate(p4_ovvvab(1:novvvab))

  novvvab = 0
  do i = 1,norb1
  do b = norb1+1,nact
     do c = norb1+1,nact
     do d = norb1+1,nact
        if (mval(ncore+i)+mval(ncore+b)==mval(ncore+c)+mval(ncore+d)) then
           novvvab = novvvab + 1
           h1_ovvvab(novvvab) = i
           p2_ovvvab(novvvab) = b
           p3_ovvvab(novvvab) = c
           p4_ovvvab(novvvab) = d
        end if
     end do
     end do
  end do
  end do

end subroutine tdcc_mkmap_ovvvab
!################################################################################
subroutine tdcc_mkmap_vvvvaa()

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : mval,nact,ncore
  use mod_cc, only : norb1,norb2
  use mod_cc2

  implicit none
  integer(c_int) :: a,b,c,d

  ! count nonzero tensors
  nvvvvaa = 0
  do a = norb1+1,nact
  do b = norb1+1,a-1
     do c = norb1+1,nact
     do d = norb1+1,c-1
        if (mval(ncore+a)+mval(ncore+b)==mval(ncore+c)+mval(ncore+d)) then
           nvvvvaa = nvvvvaa + 1
        end if
     end do
     end do
  end do
  end do
  write(6,"('tdcc_mkmap: nvvvvaa = ',2i5)") nvvvvaa,norb2**4

  allocate(p1_vvvvaa(1:nvvvvaa))
  allocate(p2_vvvvaa(1:nvvvvaa))
  allocate(p3_vvvvaa(1:nvvvvaa))
  allocate(p4_vvvvaa(1:nvvvvaa))

  nvvvvaa = 0
  do a = norb1+1,nact
  do b = norb1+1,a-1
     do c = norb1+1,nact
     do d = norb1+1,c-1
        if (mval(ncore+a)+mval(ncore+b)==mval(ncore+c)+mval(ncore+d)) then
           nvvvvaa = nvvvvaa + 1
           p1_vvvvaa(nvvvvaa) = a
           p2_vvvvaa(nvvvvaa) = b
           p3_vvvvaa(nvvvvaa) = c
           p4_vvvvaa(nvvvvaa) = d
        end if
     end do
     end do
  end do
  end do

end subroutine tdcc_mkmap_vvvvaa
!################################################################################
subroutine tdcc_mkmap_vvvvab()

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : mval,nact,ncore
  use mod_cc, only : norb1,norb2
  use mod_cc2

  implicit none
  integer(c_int) :: a,b,c,d

  ! count nonzero tensors
  nvvvvab = 0
  do a = norb1+1,nact
  do b = norb1+1,nact
     do c = norb1+1,nact
     do d = norb1+1,nact
        if (mval(ncore+a)+mval(ncore+b)==mval(ncore+c)+mval(ncore+d)) then
           nvvvvab = nvvvvab + 1
        end if
     end do
     end do
  end do
  end do
  write(6,"('tdcc_mkmap: nvvvvab = ',2i5)") nvvvvab,norb2**4

  allocate(p1_vvvvab(1:nvvvvab))
  allocate(p2_vvvvab(1:nvvvvab))
  allocate(p3_vvvvab(1:nvvvvab))
  allocate(p4_vvvvab(1:nvvvvab))

  nvvvvab = 0
  do a = norb1+1,nact
  do b = norb1+1,nact
     do c = norb1+1,nact
     do d = norb1+1,nact
        if (mval(ncore+a)+mval(ncore+b)==mval(ncore+c)+mval(ncore+d)) then
           nvvvvab = nvvvvab + 1
           p1_vvvvab(nvvvvab) = a
           p2_vvvvab(nvvvvab) = b
           p3_vvvvab(nvvvvab) = c
           p4_vvvvab(nvvvvab) = d
        end if
     end do
     end do
  end do
  end do
!  stop 'for debug @ tdcc_mkmap_vvvvab'

end subroutine tdcc_mkmap_vvvvab
!################################################################################
