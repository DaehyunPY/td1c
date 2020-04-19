!################################################################################
subroutine ormas_str()

  use, intrinsic :: iso_c_binding

  use mod_ormas, only : iprint
  use mod_ormas, only : nact, nelact, nsub
  use mod_ormas, only : nstr_alph, ndist_alph, dist_alph, nstr_alph_dist, &
       &  nstr_alph_dist_sub, arc_alph, dist_str_alph, substr_alph, onv_alph, orb_alph
  use mod_ormas, only : nstr_beta, ndist_beta, dist_beta, nstr_beta_dist, &
       &  nstr_beta_dist_sub, arc_beta, dist_str_beta, substr_beta, onv_beta, orb_beta

  implicit none
  integer(c_int) :: istr, isub

  allocate(dist_str_alph(1:2, 1:nstr_alph))
  allocate(dist_str_beta(1:2, 1:nstr_beta))
  allocate(substr_alph(1:nsub, 1:nstr_alph))
  allocate(substr_beta(1:nsub, 1:nstr_beta))
  allocate(onv_alph(1:nact, 1:nstr_alph))
  allocate(onv_beta(1:nact, 1:nstr_beta))
  allocate(orb_alph(0:nelact(1), 1:nstr_alph))
  allocate(orb_beta(0:nelact(2), 1:nstr_beta))

  if (nact == 0) return

  call ormas_str_spin(nelact(1), ndist_alph, dist_alph, nstr_alph, &
       & nstr_alph_dist, nstr_alph_dist_sub, arc_alph, dist_str_alph, &
       & substr_alph, onv_alph, orb_alph)
  call ormas_str_spin(nelact(2), ndist_beta, dist_beta, nstr_beta, &
       & nstr_beta_dist, nstr_beta_dist_sub, arc_beta, dist_str_beta, &
       & substr_beta, onv_beta, orb_beta)

  if (iprint > 0) then
     write(6, "('# ORMAS: alpha-strings', i5)") nstr_alph
     if (iprint > 1) then
        write(6, "('      #')", advance = 'no')
        write(6, "('   dist')", advance = 'no')
        write(6, "('  #dist')", advance = 'no')
        do isub = 1, nsub
           write(6, "(i7)", advance = 'no') isub
        end do
        write(6, *)
        do istr = 1, nstr_alph
           write(6, "(3i7)", advance = 'no') istr, dist_str_alph(1:2, istr)
           do isub = 1, nsub
              write(6, "(i7)", advance = 'no') substr_alph(isub, istr)
           end do
           write(6, "(2x)", advance = 'no')
           call ormas_occvec_print(6, .false., nact, onv_alph(1, istr))
           write(6, "(2x)", advance = 'no')
           if (nelact(1) > 0) call ormas_orbvec_print(.false., nelact(1), orb_alph(1, istr))
           write(6, *)
        end do
     end if

     write(6, "('# ORMAS: beta-strings', i7)") nstr_beta
     if (iprint > 1) then
        write(6, "('      #')", advance = 'no')
        write(6, "('   dist')", advance = 'no')
        write(6, "('  #dist')", advance = 'no')
        do isub = 1, nsub
           write(6, "(i7)", advance = 'no') isub
        end do
        write(6, *)
        do istr = 1, nstr_beta
           write(6, "(3i7)", advance = 'no') istr, dist_str_beta(1:2, istr)
           do isub = 1, nsub
              write(6, "(i7)", advance = 'no') substr_beta(isub, istr)
           end do
           write(6, "(2x)", advance = 'no')
           call ormas_occvec_print(6, .false., nact, onv_beta(1, istr))
           write(6, "(2x)", advance = 'no')
           if (nelact(2) > 0) call ormas_orbvec_print(.false., nelact(2), orb_beta(1, istr))
           write(6, *)
        end do
     end if
  end if

end subroutine ormas_str
!################################################################################
!################################################################################
subroutine ormas_str_spin(nel_spin, ndist_spin, dist_spin, nstr_spin, &
     & nstr_spin_dist, nstr_spin_dist_sub, arc_spin, dist_str_spin, &
     & substr_spin, onv_spin, orb_spin)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nact, nsub, norb_sub, lorb_sub

  implicit none
  integer(c_int), intent(in) :: nel_spin, ndist_spin, dist_spin(1:nsub, 1:*), &
       & nstr_spin, nstr_spin_dist(1:*), nstr_spin_dist_sub(1:nsub, 1:*), &
       & arc_spin(1:nact, 1:*)
  integer(c_int), intent(out) :: dist_str_spin(1:2, 1:*), substr_spin(1:nsub, 1:*), &
       & onv_spin(1:nact, 1:*), orb_spin(0:nel_spin, 1:*)

  integer(c_int) :: nstr_tot, istr, idist, isub, istr_dist, istr_dist_sub, rest, dim
  integer(c_int) :: nel, norb, llorb, ulorb
  integer(c_int), allocatable :: tocc(:)

  ! string --> distribution map
  nstr_tot = 0
  do idist = 1, ndist_spin
     do istr = 1, nstr_spin_dist(idist)
        nstr_tot = nstr_tot + 1
        dist_str_spin(1,nstr_tot) = idist
        dist_str_spin(2,nstr_tot) = istr
     end do
  end do
  if (nstr_tot /= nstr_spin) stop 'error in det_dist_str.'

  ! string --> substring map
  do istr = 1, nstr_spin
     idist = dist_str_spin(1, istr)
     istr_dist = dist_str_spin(2, istr)

     rest = istr_dist - 1
     dim = nstr_spin_dist(idist)
     do isub = 1, nsub
        dim = dim / nstr_spin_dist_sub(isub, idist)
        do istr_dist_sub = 1, nstr_spin_dist_sub(isub, idist)
           if (rest < dim * istr_dist_sub) then
              substr_spin(isub, istr) = istr_dist_sub
              rest = rest - dim * (istr_dist_sub - 1)
              exit
           end if
        end do
     end do
  end do

  ! occupation number vector
  allocate(tocc(1:nact))
  do istr = 1, nstr_spin
     idist = dist_str_spin(1, istr)
     do isub = 1, nsub
        nel = dist_spin(isub, idist)
        norb = norb_sub(isub)
        istr_dist_sub = substr_spin(isub, istr)
        call ormas_occvec(nel, norb, istr_dist_sub, arc_spin, tocc)

        llorb = lorb_sub(1, isub)
        ulorb = lorb_sub(2, isub)
        onv_spin(llorb:ulorb, istr) = tocc(1:norb)
     end do
  end do
  deallocate(tocc)

  ! orbital index vector
  do istr = 1, nstr_spin
     call ormas_orbvec(nact, onv_spin(1, istr), orb_spin(0, istr))
  end do

end subroutine ormas_str_spin
!################################################################################
