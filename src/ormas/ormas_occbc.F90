!################################################################################
subroutine ormas_occbc()

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : iprint, nact, nsub, norb_sub, lorb_sub, &
       & sub_orb, min_sub, max_sub, min_sub_alph, max_sub_alph, &
       & min_sub_beta, max_sub_beta

  implicit none
  integer(c_int) :: isub, iact

  allocate(sub_orb(1:nact))
  do isub = 1, nsub
     do iact = lorb_sub(1, isub), lorb_sub(2, isub)
        sub_orb(iact) = isub
     end do
  end do
  call ormas_occbcx()

  if (iprint > 0) then
     write(6, "('# ORMAS: orbital subspaces')")
     write(6, "('#                 norb     llorb     ulorb')")
     do isub = 1, nsub
        write(6, "(4i10)") isub, norb_sub(isub), lorb_sub(1:2, isub)
     end do

     write(6, "('# ORMAS: orb --> sub map')")
     write(6, "('#                 orb       sub')")
     do iact = 1, nact
        isub = sub_orb(iact)
        write(6, "(10x, 2i10)") iact, isub
     end do

     write(6, "('# ORMAS: occupation number boundaries')")
     write(6, "('#                 min       max      min_a     max_a     min_b     max_b')")
     do isub = 1, nsub
        write(6, "(7i10)") isub, &
             & min_sub(isub), max_sub(isub), &
             & min_sub_alph(isub), max_sub_alph(isub), &
             & min_sub_beta(isub), max_sub_beta(isub)
     end do
  end if

end subroutine ormas_occbc
!################################################################################
!################################################################################
subroutine ormas_occbcx()

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nelact, nsub, norb_sub, min_sub, max_sub
  use mod_ormas, only : min_sub_alph, max_sub_alph, min_sub_beta, max_sub_beta

  implicit none
  integer(c_int) :: isub, jsub, rmax

  do isub = 1, nsub
     max_sub_alph(isub) = min(norb_sub(isub), min(max_sub(isub), nelact(1)))
     max_sub_beta(isub) = min(norb_sub(isub), min(max_sub(isub), nelact(2)))
  end do

  do isub = 1, nsub
     rmax = nelact(1)
     do jsub = 1, nsub
        if (jsub /= isub) then
           rmax = rmax - max_sub_alph(jsub)
        end if
     end do
     min_sub_alph(isub) = max(0, max(min_sub(isub) - max_sub_beta(isub), rmax))

     rmax = nelact(2)
     do jsub = 1, nsub
        if (jsub /= isub) then
           rmax = rmax - max_sub_beta(jsub)
        end if
     end do
     min_sub_beta(isub) = max(0, max(min_sub(isub) - max_sub_alph(isub), rmax))
  end do  

end subroutine ormas_occbcx
!################################################################################
