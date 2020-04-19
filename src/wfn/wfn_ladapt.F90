!///////////////////////////////////////////////////////////////////////
subroutine wfn_ladapt(orb, cic)

  use, intrinsic :: iso_c_binding
  implicit none
  complex(c_double_complex), intent(inout) :: orb(1:*)
  complex(c_double_complex), intent(inout) :: cic(1:*)

  call wfn_ladapto(orb)
  call wfn_ladaptc(cic)

end subroutine wfn_ladapt
!///////////////////////////////////////////////////////////////////////
subroutine wfn_ladapto(orb)

  use, intrinsic :: iso_c_binding

  use mod_rad, only : nrad
  use mod_sph, only : lmax1
  use mod_const, only : czero
  use mod_bas, only : lval
  use mod_ormas, only : nfun

  implicit none
  complex(c_double_complex), intent(inout) :: orb(1:(nrad-1), 0:lmax1, 1:nfun)
  integer(c_int) :: l, irad, ifun

  do ifun = 1, nfun
     do l = 0, lmax1
        if (l.ne.lval(ifun)) then
!       if (mod(l-lval(ifun),2).ne.0) then
           do irad = 1, nrad - 1
              orb(irad, l, ifun) = czero
           end do
        end if
     end do
  end do

end subroutine wfn_ladapto
!///////////////////////////////////////////////////////////////////////
subroutine wfn_ladaptc(cic)

  use, intrinsic :: iso_c_binding

  use mod_rad, only : nrad
  use mod_sph, only : lmax1
  use mod_const, only : czero
  use mod_bas, only : lval, ltot
  use mod_ormas, only : nfun,nelact,ncore,nact
  use mod_ormas, only : orb_alph,orb_beta,lcic,ndetx,tdcc
  use mod_ormas, only : mapr_detx
  !nstr_alph,nstr_beta,llstr_alph_beta,nstr_alph_beta

  implicit none
  complex(c_double_complex), intent(inout) :: cic(1:lcic)
  integer(c_int) :: idet, istr, jstr, suma, sumb, iela, ielb

  if (nact == 0) return
  if (ltot /= 0) return
  !DEBUG
  !write(6,"('wfn_ladaptc is disabled.')")
  return
  !DEBUG

  do idet = 1, ndetx
     jstr = mapr_detx(1,idet)
     istr = mapr_detx(2,idet)
     suma = 0
     do iela = 1, nelact(1)
        suma = suma + lval(ncore+orb_alph(iela,jstr))
     end do
     sumb = 0
     do ielb = 1, nelact(2)
        sumb = sumb + lval(ncore+orb_beta(ielb,istr))
     end do
     if (mod(suma+sumb,2) .ne. 0) then
        cic(idet) = czero
        if (tdcc) cic(ndetx+idet) = czero
     end if
  end do

!old  do istr = 1, nstr_beta
!old     sumb = 0
!old     do ielb = 1, nelact(2)
!old        sumb = sumb + lval(ncore+orb_beta(ielb,istr))
!old     end do
!old     do jstr = llstr_alph_beta(istr), llstr_alph_beta(istr)+nstr_alph_beta(istr)-1
!old        suma = 0
!old        do iela = 1, nelact(1)
!old           suma = suma + lval(ncore+orb_alph(iela,jstr))
!old        end do
!old        if (mod(suma+sumb,2) .ne. 0) then
!old           cic(mapf_detx(jstr,istr)) = czero
!old        end if
!old     end do
!old  end do

end subroutine wfn_ladaptc
!///////////////////////////////////////////////////////////////////////
subroutine wfn_ladapt_core(orb)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad
  use mod_sph, only : lmax1
  use mod_bas, only : lval
  use mod_ormas, only : ncore, nfun
  use mod_const, only : czero

  implicit none
  complex(c_double_complex), intent(inout) :: orb(1:(nrad-1), 0:lmax1, 1:nfun)

  integer(c_int) :: ifun, l, li, irad

  do ifun = 1, ncore
     li = lval(ifun)
     do l = 0, lmax1
        if (l .ne. li) then
           do irad = 1, nrad - 1
              orb(irad, l, ifun) = czero
           end do
        end if
     end do
  end do

end subroutine wfn_ladapt_core
!///////////////////////////////////////////////////////////////////////
