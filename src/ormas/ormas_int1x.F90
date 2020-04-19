!################################################################################
subroutine ormas_int1x()

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : iprint
  use mod_ormas, only : nelact, nact
  use mod_ormas, only : nstr_alph, ndist_alph, min_sub_alph, max_sub_alph, &
       & arc_alph, onv_alph, dist_alph, nstr_alph_dist, nstr_alph_dist_sub, dist_str_alph, &
       & n1x_alph, p1x_alph, h1x_alph, eq1x_alph, sgn1x_alph
  use mod_ormas, only : nstr_beta, ndist_beta, min_sub_beta, max_sub_beta, &
       & arc_beta, onv_beta, dist_beta, nstr_beta_dist, nstr_beta_dist_sub, dist_str_beta, &
       & n1x_beta, p1x_beta, h1x_beta, eq1x_beta, sgn1x_beta
  use mod_ormas, only : n1x_m_alph, map1x_m_alph, n1x_m_beta, map1x_m_beta

  implicit none
  integer(c_long) :: na2

  na2 = nact * nact
  allocate(n1x_alph(-3:1, 1:nstr_alph))
  allocate(p1x_alph(1:na2, 1:nstr_alph))
  allocate(h1x_alph(1:na2, 1:nstr_alph))
  allocate(eq1x_alph(1:na2, 1:nstr_alph))
  allocate(sgn1x_alph(1:na2, 1:nstr_alph))

  allocate(n1x_beta(-3:1, 1:nstr_beta))
  allocate(p1x_beta(1:na2, 1:nstr_beta))
  allocate(h1x_beta(1:na2, 1:nstr_beta))
  allocate(eq1x_beta(1:na2, 1:nstr_beta))
  allocate(sgn1x_beta(1:na2, 1:nstr_beta))

  call ormas_int1x_spin(nelact(1), min_sub_alph, max_sub_alph, ndist_alph, dist_alph, &
       & nstr_alph, nstr_alph_dist, nstr_alph_dist_sub, dist_str_alph, arc_alph, onv_alph, n1x_alph, &
       & p1x_alph, h1x_alph, eq1x_alph, sgn1x_alph)
  call ormas_int1x_spin(nelact(2), min_sub_beta, max_sub_beta, ndist_beta, dist_beta, &
       & nstr_beta, nstr_beta_dist, nstr_beta_dist_sub, dist_str_beta, arc_beta, onv_beta, n1x_beta, &
       & p1x_beta, h1x_beta, eq1x_beta, sgn1x_beta)

  if (iprint > 2) then
     call ormas_int1x_print()
  end if

end subroutine ormas_int1x
!################################################################################
!################################################################################
subroutine ormas_int1x_spin(nel, min_sub, max_sub, ndist, dist, &
     & nstr, nstr_dist, nstr_dist_sub, dist_str, arc, onv, n1x, p1x, h1x, eq1x, sgn1x)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : iprint
  use mod_ormas, only : ncore, nact, nsub, lorb_sub

  implicit none
  integer(c_long), intent(in) :: nel, nstr, ndist
  integer(c_long), intent(in) :: min_sub(1:*)
  integer(c_long), intent(in) :: max_sub(1:*)
  integer(c_long), intent(in) :: arc(1:nact, 1:*)
  integer(c_long), intent(in) :: onv(1:nact, 1:*)
  integer(c_long), intent(in) :: dist(1:nsub, 1:*)
  integer(c_long), intent(in) :: nstr_dist(1:*)
  integer(c_long), intent(in) :: nstr_dist_sub(1:nsub, 1:*)
  integer(c_long), intent(in) :: dist_str(1:2, 1:*)
  integer(c_long), intent(out) :: n1x(-3:1, 1:*)
  integer(c_long), intent(out) :: p1x(1:nact*nact, 1:*)
  integer(c_long), intent(out) :: h1x(1:nact*nact, 1:*)
  integer(c_long), intent(out) :: eq1x(1:nact*nact, 1:*)
  integer(c_long), intent(out) :: sgn1x(1:nact*nact, 1:*)

  integer(c_long), external :: ormas_map1x, ormas_sgn1x
  integer(c_long), allocatable :: tocc(:)
  integer(c_long) :: istr, ifun, jfun, idist
  integer(c_long) :: isub, jsub, neli, nelj, lli, uli, llj, ulj, tmp

  allocate(tocc(1:nact))

  do istr = 1, nstr
     n1x(-3:1, istr) = 0
     idist = dist_str(1, istr)
     do jsub = 1, nsub ! hole
        nelj = dist(jsub, idist)
        llj = lorb_sub(1, jsub)
        ulj = lorb_sub(2, jsub)

        do isub = 1, nsub ! particle
           neli = dist(isub, idist)
           lli = lorb_sub(1, isub)
           uli = lorb_sub(2, isub)

           if (jsub == isub) then
              tmp = 0
           else if (nelj == min_sub(jsub) .and. neli == max_sub(isub)) then
              tmp = -3
           else if (neli == max_sub(isub)) then
              tmp = -2
           else if (nelj == min_sub(jsub)) then
              tmp = -1
           else
              tmp = 0
           end if

           do jfun = llj, ulj ! hole
              if(onv(jfun, istr) == 0) cycle

              do ifun = lli, uli ! particle
                 if(ifun /= jfun .and. onv(ifun, istr) /= 0) cycle
                 !M-adapt
                 !if(mval(ncore+ifun) /= mval(ncore+jfun)) cycle
                 !M-adapt

                 tocc(1:nact) = onv(1:nact, istr)
                 tocc(ifun) = 1
                 if (ifun /= jfun) tocc(jfun) = 0

                 n1x(1,   istr) = n1x(1,   istr) + 1
                 n1x(tmp, istr) = n1x(tmp, istr) + 1

                 p1x(n1x(1, istr), istr) = ifun
                 h1x(n1x(1, istr), istr) = jfun
                 if (tmp == 0) then
                    eq1x (n1x(1, istr), istr) = ormas_map1x(nel, nact, ndist, dist, &
                         &  nstr_dist, nstr_dist_sub, arc, tocc)
                    sgn1x(n1x(1, istr), istr) = ormas_sgn1x(nact, ifun, jfun, onv(1,istr))
                 else
                    eq1x (n1x(1, istr), istr) = tmp
                    sgn1x(n1x(1, istr), istr) = ormas_sgn1x(nact, ifun, jfun, onv(1,istr))
                 end if
              end do
           end do
        end do
     end do

     call ormas_int1x_order(n1x(-3, istr), p1x(1, istr), h1x(1, istr), &
          &  eq1x(1, istr), sgn1x(1, istr))
  end do

  deallocate(tocc)

end subroutine ormas_int1x_spin
!################################################################################
!################################################################################
subroutine ormas_int1x_order(n1x, p1x, h1x, eq1x, sgn1x)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nact

  implicit none
  integer(c_long), intent(in) :: n1x(-3:1)
  integer(c_long), intent(inout) :: p1x(1:nact*nact)
  integer(c_long), intent(inout) :: h1x(1:nact*nact)
  integer(c_long), intent(inout) :: eq1x(1:nact*nact)
  integer(c_long), intent(inout) :: sgn1x(1:nact*nact)

  integer(c_long) :: offs(-3:0), n1x_type(-3:0), na2, i1x, type
  integer(c_long), allocatable :: p1x_tmp(:), h1x_tmp(:), eq1x_tmp(:), sgn1x_tmp(:)

  offs(0) = 0
  offs(-1) = offs(0) + n1x(0)
  offs(-2) = offs(-1) + n1x(-1)
  offs(-3) = offs(-2) + n1x(-2)

  na2 = nact * nact
  allocate(p1x_tmp(1:na2))
  allocate(h1x_tmp(1:na2))
  allocate(eq1x_tmp(1:na2))
  allocate(sgn1x_tmp(1:na2))

  p1x_tmp(1:na2) = p1x(1:na2)
  h1x_tmp(1:na2) = h1x(1:na2)
  eq1x_tmp(1:na2) = eq1x(1:na2)
  sgn1x_tmp(1:na2) = sgn1x(1:na2)

  n1x_type(-3:0) = 0
  do i1x = 1, n1x(1)
     type = min(0, eq1x_tmp(i1x))
     n1x_type(type) = n1x_type(type) + 1
     p1x(offs(type) + n1x_type(type)) = p1x_tmp(i1x)
     h1x(offs(type) + n1x_type(type)) = h1x_tmp(i1x)
     eq1x(offs(type) + n1x_type(type)) = eq1x_tmp(i1x)
     sgn1x(offs(type) + n1x_type(type)) = sgn1x_tmp(i1x)
  end do

  deallocate(sgn1x_tmp)
  deallocate(eq1x_tmp)
  deallocate(h1x_tmp)
  deallocate(p1x_tmp)

end subroutine ormas_int1x_order
!################################################################################
!################################################################################
subroutine ormas_int1x_print()

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : iprint
  use mod_ormas, only : nact
  use mod_ormas, only : nstr_alph, onv_alph, n1x_alph, p1x_alph, h1x_alph, &
       & eq1x_alph, sgn1x_alph
  use mod_ormas, only : nstr_beta, onv_beta, n1x_beta, p1x_beta, h1x_beta, &
       & eq1x_beta, sgn1x_beta

  implicit none

  if (iprint > 2) then
     write(6, "(' # ORMAS: alpha int1x')")
     call ormas_int1x_print_spin(nact, nstr_alph, onv_alph, n1x_alph, &
          & p1x_alph, h1x_alph, eq1x_alph, sgn1x_alph)
     write(6, "(' # ORMAS: beta int1x')")
     call ormas_int1x_print_spin(nact, nstr_beta, onv_beta, n1x_beta, &
          & p1x_beta, h1x_beta, eq1x_beta, sgn1x_beta)
  end if

end subroutine ormas_int1x_print
!################################################################################
!################################################################################
subroutine ormas_int1x_print_spin(nact, nstr, onv, n1x, p1x, h1x, &
     & eq1x, sgn1x)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : iprint

  implicit none
  integer(c_long), intent(in) :: nact, nstr
  integer(c_long), intent(in) :: onv(1:nact, 1:*)
  integer(c_long), intent(in) :: n1x(-3:1, 1:*)
  integer(c_long), intent(in) :: p1x(1:nact*nact, 1:*)
  integer(c_long), intent(in) :: h1x(1:nact*nact, 1:*)
  integer(c_long), intent(in) :: eq1x(1:nact*nact, 1:*)
  integer(c_long), intent(in) :: sgn1x(1:nact*nact, 1:*)

  integer(c_long), allocatable :: tocc(:)
  integer(c_long) :: istr, i1x, ifun, jfun

  allocate(tocc(1:nact))

  do istr = 1, nstr
     do i1x = 1, n1x(1, istr)
        ifun = p1x(i1x, istr) ! particle
        jfun = h1x(i1x, istr) ! hole

        tocc(1:nact) = onv(1:nact, istr)
        tocc(ifun) = 1
        if (ifun /= jfun) tocc(jfun) = 0

        write(6,"(2i5)",advance='no') istr, i1x
        write(6,"(2x)",advance='no')

        call ormas_occvec_print(6, .false., nact, onv(1,istr))
        write(6,"(' --(')",advance='no')
        write(6,"(2i5)",advance='no') h1x(i1x,istr), p1x(i1x,istr)
        write(6,"(')-> ')",advance='no')
        call ormas_occvec_print(6, .false., nact, tocc)
        write(6,"(2x)",advance='no')
!                    write(6,"(';',2x)",advance='no')
!                    call det_string(nact, onv(1,istr), tstr)
!                    call det_print_string(.false., nel, tstr)
!                    write(6,"(' -->',x)",advance='no')
!                    call det_string(nact, tocc, tstr)
!                    call det_print_string(.false., nel, tstr)
        write(6,"(i5)", advance='no') eq1x(i1x, istr)
        write(6,"(2x)", advance='no')
        write(6,"(i5)", advance='no') sgn1x(i1x, istr)
        write(6,*)
     end do
  end do

  deallocate(tocc)

end subroutine ormas_int1x_print_spin
!################################################################################
