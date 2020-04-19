!################################################################################
subroutine ormas_int1xr()

!OLD  use, intrinsic :: iso_c_binding
!OLD  use mod_ormas, only : iprint
!OLD  use mod_ormas, only : nelact, nact
!OLD  use mod_ormas, only : nstr_alph, ndist_alph, min_sub_alph, max_sub_alph, &
!OLD       & arc_alph, onv_alph, dist_alph, nstr_alph_dist, nstr_alph_dist_sub, dist_str_alph, &
!OLD       & n1xr_alph, r1xr_alph, l1xr_alph, sgn1xr_alph
!OLD  use mod_ormas, only : nstr_beta, ndist_beta, min_sub_beta, max_sub_beta, &
!OLD       & arc_beta, onv_beta, dist_beta, nstr_beta_dist, nstr_beta_dist_sub, dist_str_beta, &
!OLD       & n1xr_beta, r1xr_beta, l1xr_beta, sgn1xr_beta

  implicit none
  stop "ormas_int1xr no longer supported."

!OLD  allocate(n1xr_alph  (             1:nact, 1:nact))
!OLD  allocate(r1xr_alph  (1:nstr_alph, 1:nact, 1:nact))
!OLD  allocate(l1xr_alph  (1:nstr_alph, 1:nact, 1:nact))
!OLD  allocate(sgn1xr_alph(1:nstr_alph, 1:nact, 1:nact))
!OLD  allocate(n1xr_beta  (             1:nact, 1:nact))
!OLD  allocate(r1xr_beta  (1:nstr_beta, 1:nact, 1:nact))
!OLD  allocate(l1xr_beta  (1:nstr_beta, 1:nact, 1:nact))
!OLD  allocate(sgn1xr_beta(1:nstr_beta, 1:nact, 1:nact))
!OLD
!OLD  call ormas_int1xr_spin(nelact(1), min_sub_alph, max_sub_alph, ndist_alph, dist_alph, &
!OLD       & nstr_alph, nstr_alph_dist, nstr_alph_dist_sub, dist_str_alph, arc_alph, onv_alph, n1xr_alph, &
!OLD       & r1xr_alph, l1xr_alph, sgn1xr_alph)
!OLD  call ormas_int1xr_spin(nelact(2), min_sub_beta, max_sub_beta, ndist_beta, dist_beta, &
!OLD       & nstr_beta, nstr_beta_dist, nstr_beta_dist_sub, dist_str_beta, arc_beta, onv_beta, n1xr_beta, &
!OLD       & r1xr_beta, l1xr_beta, sgn1xr_beta)
!OLD
!OLD  if (iprint > 2) then
!OLD     call ormas_int1xr_print()
!OLD  end if

end subroutine ormas_int1xr
!################################################################################
!################################################################################
subroutine ormas_int1xr_spin(nel, min_sub, max_sub, ndist, dist, &
     & nstr, nstr_dist, nstr_dist_sub, dist_str, arc, onv, n1xr, r1xr, l1xr, sgn1xr)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : iprint
  use mod_ormas, only : ncore, nact, nsub, lorb_sub

  implicit none
  integer(c_int), intent(in) :: nel, nstr, ndist
  integer(c_int), intent(in) :: min_sub(1:*)
  integer(c_int), intent(in) :: max_sub(1:*)
  integer(c_int), intent(in) :: arc(1:nact, 1:*)
  integer(c_int), intent(in) :: onv(1:nact, 1:*)
  integer(c_int), intent(in) :: dist(1:nsub, 1:*)
  integer(c_int), intent(in) :: nstr_dist(1:*)
  integer(c_int), intent(in) :: nstr_dist_sub(1:nsub, 1:*)
  integer(c_int), intent(in) :: dist_str(1:2, 1:*)
  integer(c_int), intent(out) :: n1xr(1:nact, 1:nact)
  integer(c_int), intent(out) :: r1xr(1:nstr, 1:nact, 1:nact)
  integer(c_int), intent(out) :: l1xr(1:nstr, 1:nact, 1:nact)
  integer(c_int), intent(out) :: sgn1xr(1:nstr, 1:nact, 1:nact)

  integer(c_int), external :: ormas_map1x, ormas_sgn1x
  integer(c_int), allocatable :: tocc(:)
  integer(c_int) :: istr, ifun, jfun, idist
  integer(c_int) :: isub, jsub, neli, nelj, lli, uli, llj, ulj

  allocate(tocc(1:nact))

  do isub = 1, nsub ! particle
     lli = lorb_sub(1, isub)
     uli = lorb_sub(2, isub)
     do ifun = lli, uli
        do jsub = 1, nsub ! hole
           llj = lorb_sub(1, jsub)
           ulj = lorb_sub(2, jsub)
           do jfun = llj, ulj

              n1xr(jfun, ifun) = 0
              do istr = 1, nstr
                 idist = dist_str(1, istr)
                 neli = dist(isub, idist)
                 nelj = dist(jsub, idist)

                 if (jsub /= isub .and. (nelj == min_sub(jsub) .or. neli == max_sub(isub))) cycle
                 if(onv(jfun, istr) == 0 .or. (ifun /= jfun .and. onv(ifun, istr) /= 0)) cycle
                 !M-adapt
                 !if(mval(ncore+ifun) /= mval(ncore+jfun)) cycle
                 !M-adapt

                 tocc(1:nact) = onv(1:nact, istr)
                 tocc(ifun) = 1                   ! particle
                 if (ifun /= jfun) tocc(jfun) = 0 ! hole

                 n1xr(jfun, ifun) = n1xr(jfun, ifun) + 1
                 r1xr(n1xr(jfun, ifun), jfun, ifun) = istr

                 l1xr(n1xr(jfun, ifun), jfun, ifun) = ormas_map1x(nel, nact, ndist, dist, &
                      & nstr_dist, nstr_dist_sub, arc, tocc)
                 sgn1xr(n1xr(jfun, ifun), jfun, ifun) = ormas_sgn1x(nact, ifun, &
                      & jfun, onv(1,istr))

              end do
           end do
        end do
     end do
  end do

  deallocate(tocc)

end subroutine ormas_int1xr_spin
!################################################################################
!################################################################################
subroutine ormas_int1xr_print()

  use, intrinsic :: iso_c_binding
!OLD  use mod_ormas, only : iprint
!OLD  use mod_ormas, only : nelact, nact
!OLD  use mod_ormas, only : nstr_alph, onv_alph, n1xr_alph, r1xr_alph, l1xr_alph, &
!OLD       & sgn1xr_alph
!OLD  use mod_ormas, only : nstr_beta, onv_beta, n1xr_beta, r1xr_beta, l1xr_beta, &
!OLD       & sgn1xr_beta

  implicit none
  stop "ormas_int1xr_print no longer supported."

!OLD  if (iprint > 2) then
!OLD     write(6, "(' # ORMAS: alpha int1xr')")
!OLD     call ormas_int1xr_print_spin(nelact(1), nact, nstr_alph, onv_alph, &
!OLD          & n1xr_alph, r1xr_alph, l1xr_alph, sgn1xr_alph)
!OLD     write(6, "(' # ORMAS: beta int1xr')")
!OLD     call ormas_int1xr_print_spin(nelact(1), nact, nstr_alph, onv_alph, &
!OLD          & n1xr_alph, r1xr_alph, l1xr_alph, sgn1xr_alph)
!OLD!  call ormas_int1xr_print_spin(nelact(2), nact, nstr_beta, onv_beta, &
!OLD!       & n1xr_beta, r1xr_beta, l1xr_beta, sgn1xr_beta)
!OLD  end if

end subroutine ormas_int1xr_print
!################################################################################
!################################################################################
subroutine ormas_int1xr_print_spin(nel, nact, nstr, onv, n1xr, r1xr, l1xr, &
     & sgn1xr)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : iprint

  implicit none
  integer(c_int), intent(in) :: nel, nact, nstr
  integer(c_int), intent(in) :: onv(1:nact, *)
  integer(c_int), intent(in) :: n1xr  (        1:nact, 1:nact)
  integer(c_int), intent(in) :: r1xr  (1:nstr, 1:nact, 1:nact)
  integer(c_int), intent(in) :: l1xr  (1:nstr, 1:nact, 1:nact)
  integer(c_int), intent(in) :: sgn1xr(1:nstr, 1:nact, 1:nact)

  integer(c_int), allocatable :: tocc(:)
  integer(c_int) :: istr, ifun, jfun, i1xr

  allocate(tocc(1:nact))

  do ifun = 1, nact
     do jfun = 1, nact
        do i1xr = 1, n1xr(jfun, ifun)
           istr = r1xr(i1xr, jfun, ifun)
           
           tocc(1:nact) = onv(1:nact, istr)
           tocc(ifun) = 1                   ! particle
           if (ifun /= jfun) tocc(jfun) = 0 ! hole
           
           write(6,"(3i5)",advance='no') jfun, ifun, i1xr
           write(6,"(2x)",advance='no')
           
           write(6,"(i5)",advance='no') istr
           write(6,"('(')",advance='no')
           call ormas_occvec_print(6, .false., nact, onv(1,istr))
           write(6,"(')')",advance='no')
           write(6,"(' --> ')",advance='no')
           write(6,"(i5)",advance='no') l1xr(i1xr, jfun, ifun)
           write(6,"('(')",advance='no')
           call ormas_occvec_print(6, .false., nact, tocc)
           write(6,"(')')",advance='no')
           write(6,"(2x)",advance='no')
           write(6,"(i5)") sgn1xr(i1xr, jfun, ifun)
        end do
     end do
  end do

!  do istr = 1, nstr
!     n1x(istr) = 0
!
!     do jfun = 1, nact ! anihilation
!        if(onv(jfun, istr)) then
!
!           do ifun = 1, nact ! creation
!              if(ifun==jfun .or. .not.onv(ifun, istr)) then
!
!                 n1x(istr) = n1x(istr) + 1
!                 p1x(n1x(istr), istr) = ifun ! particle
!                 h1x(n1x(istr), istr) = jfun ! hole
!
!                 tocc(1:nact) = onv(1:nact, istr)
!                 tocc(ifun) = .true.                    ! particle
!                 if (ifun /= jfun) tocc(jfun) = .false. ! hole
!
!                 eq1x (n1x(istr), istr) = det_map1(nact, tocc, wgt(1,1))
!                 sgn1x(n1x(istr), istr) = det_sgn1x(nact, ifun, jfun, onv(1,istr))
!
!                 if (iprint > 2) then
!                    write(6,"(2i5)",advance='no') istr, n1x(istr)
!                    write(6,"(2x)",advance='no')
!
!                    call det_print_occ(.false., nact, onv(1,istr))
!                    write(6,"(' --(')",advance='no')
!                    write(6,"(2i5)",advance='no') h1x(n1x(istr),istr), p1x(n1x(istr),istr)
!                    write(6,"(')-> ')",advance='no')
!                    call det_print_occ(.false., nact, tocc)
!                    write(6,"(2x)",advance='no')
!!                    write(6,"(';',2x)",advance='no')
!!                    call det_string(nact, onv(1,istr), tstr)
!!                    call det_print_string(.false., nel, tstr)
!!                    write(6,"(' -->',x)",advance='no')
!!                    call det_string(nact, tocc, tstr)
!!                    call det_print_string(.false., nel, tstr)
!                    write(6,"(i5)",advance='no') eq1x(n1x(istr), istr)
!                    write(6,"(2x)",advance='no')
!                    write(6,"(i5)") sgn1x(n1x(istr), istr)
!                 end if
!              end if
!           end do
!        end if
!     end do
!  end do

  deallocate(tocc)

end subroutine ormas_int1xr_print_spin
!################################################################################
