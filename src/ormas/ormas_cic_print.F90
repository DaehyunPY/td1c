!################################################################################
subroutine ormas_cic_print(cic, fname)

  use, intrinsic :: iso_c_binding

  use mod_bas, only : mval
  use mod_ormas, only : nact, ncore, nelact, ntot_alph_beta
  use mod_ormas, only : nstr_alph, dist_str_alph, onv_alph, orb_alph
  use mod_ormas, only : nstr_beta, dist_str_beta, onv_beta, orb_beta
  use mod_ormas, only : mval_alph, mval_beta, tdcc, lcic, ndetx
  use mod_ormas, only : llstr_alph_beta, nstr_alph_beta

  implicit none
  character(*), intent(in) :: fname
  complex(c_double_complex), intent(in) :: cic(1:lcic)
  integer(c_int) :: istr, jstr, ifun, jfun
  integer(c_int), parameter :: ioout = 1

  if (nact == 0) return
  open(unit=ioout, file=trim(fname), status='unknown', form='formatted')
  call ormas_cic_printx(cic, ioout)
  close(unit=ioout)

  !debug
  if (tdcc)  call tdcc_print(cic)
  !debug

end subroutine ormas_cic_print
!################################################################################
subroutine ormas_cic_printx(cic, ioout)

  use, intrinsic :: iso_c_binding

  use mod_bas, only : mval
  use mod_ormas, only : nact, ncore, nelact, ntot_alph_beta
  use mod_ormas, only : nstr_alph, dist_str_alph, onv_alph, orb_alph
  use mod_ormas, only : nstr_beta, dist_str_beta, onv_beta, orb_beta
  use mod_ormas, only : mval_alph, mval_beta, tdcc, lcic, ndetx
  use mod_ormas, only : llstr_alph_beta, nstr_alph_beta

  implicit none
  integer(c_int), intent(in) :: ioout
  complex(c_double_complex), intent(in) :: cic(1:lcic)
  integer(c_int) :: istr, jstr, ifun, jfun

  if (nact == 0) return
  do istr = 1, nstr_beta
!     do jstr = 1, nstr_alph
     do jstr = llstr_alph_beta(istr), llstr_alph_beta(istr)+nstr_alph_beta(istr)-1
        write(ioout, "(3i7)", advance = 'no') &
             dist_str_alph(1,jstr)-1,dist_str_beta(1,istr)-1, &
             dist_str_alph(1,jstr) + dist_str_beta(1,istr)-2
        write(ioout, "(3i7)", advance = 'no') mval_alph(jstr),mval_beta(istr),mval_alph(jstr)+mval_beta(istr)
        write(ioout, "(5x)", advance = 'no')
        call ormas_occvec_print(ioout, .false., nact,onv_alph(1,jstr))
        write(ioout, "(' x ')", advance = 'no')
        call ormas_occvec_print(ioout, .false., nact,onv_beta(1,istr))
        write(ioout, "(2f20.10)", advance = 'no') cic(ntot_alph_beta(istr)+jstr)
        if (tdcc) write(ioout, "(2f20.10)", advance = 'no') cic(ndetx+ntot_alph_beta(istr)+jstr)
!debug
!        write(6,"(': dist = ',3i3)") dist_alph(isub,jdist),dist_beta(isub,idist),dist_alph(isub,jdist)+dist_beta(isub,idist)
!debug
        write(ioout,*)
     end do
  end do

  !debug
  !if (tdcc)  call tdcc_print(cic)
  !debug

end subroutine ormas_cic_printx
!################################################################################
subroutine ormas_cic_print_old(cic, fname)

  use, intrinsic :: iso_c_binding

  use mod_const, only : czero
  use mod_ormas, only : nact, ncore, nelact, ntot_alph_beta, det_allowed
  use mod_ormas, only : nstr_alph, dist_str_alph, onv_alph, orb_alph
  use mod_ormas, only : nstr_beta, dist_str_beta, onv_beta, orb_beta
  use mod_ormas, only : llstr_alph_beta, nstr_alph_beta, map1to2_alph, map1to2_beta

  implicit none
  character(*), intent(in) :: fname
  complex(c_double_complex), intent(in) :: cic(1:*)
  integer(c_int) :: istr1, jstr1, istr2, jstr2, ifun, jfun

  if (nact == 0) return
  stop "ormas_cic_print_old no longer valid"
!old  open(unit=1, file=trim(fname), status='unknown', form='formatted')
!old  do istr1 = 1, nstr_beta
!old     istr2 = map1to2_beta(istr1)
!old     do jstr1 = 1, nstr_alph
!old        jstr2 = map1to2_alph(jstr1)
!old        if (det_allowed(dist_str_alph(1,jstr2),dist_str_beta(1,istr2)) == 0) cycle
!old        write(1, "(3i7)", advance = 'no') &
!old             dist_str_alph(1,jstr2)-1,dist_str_beta(1,istr2)-1, &
!old             dist_str_alph(1,jstr2) + dist_str_beta(1,istr2)-2
!old        write(1, "(5x)", advance = 'no')
!old        call ormas_occvec_print(1, .false., nact,onv_alph(1,jstr2))
!old        write(1, "(' x ')", advance = 'no')
!old        call ormas_occvec_print(1, .false., nact,onv_beta(1,istr2))
!old        if (mapf_detx(jstr2,istr2) .ne. 0) then
!old           write(1, "(2f20.10)") cic(mapf_detx(jstr2,istr2))
!old        else
!old           write(1, "(2f20.10)") czero
!old        end if
!old     end do
!old  end do
!old  close(unit=1)

end subroutine ormas_cic_print_old
!################################################################################
