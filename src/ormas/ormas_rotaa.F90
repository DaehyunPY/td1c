!################################################################################
subroutine ormas_rotaa()

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : iprint
  use mod_ormas, only : ncore, nact, nfun, mval, nsub, norb_sub, lorb_sub, nrotaa, rotaa_mapf, rotaa_mapb

  implicit none
  integer(c_long) :: isub, jsub, iact, jact, irot, npair, check

  ! number of inter-subactive space rotatins
  nrotaa = 0
  !do isub = 1, nsub
  !   do jsub = 1, isub - 1
  !      nrotaa = nrotaa + norb_sub(jsub) * norb_sub(isub)
  !   end do
  !end do
  do isub = 1, nsub
     do jsub = 1, isub - 1
        do iact = lorb_sub(1, isub), lorb_sub(2, isub)
           do jact = lorb_sub(1, jsub), lorb_sub(2, jsub)
              if (mval(ncore+iact) == mval(ncore+jact)) nrotaa = nrotaa + 1
           end do
        end do
     end do
  end do

  allocate(rotaa_mapf(1:nact, 1:nact))
  allocate(rotaa_mapb(1:nrotaa, 1:2))

  npair = 0
  do isub = 1, nsub
     do jsub = 1, isub - 1
        do iact = lorb_sub(1, isub), lorb_sub(2, isub)
           do jact = lorb_sub(1, jsub), lorb_sub(2, jsub)
              if (mval(ncore+iact) == mval(ncore+jact)) then
                 npair = npair + 1
                 rotaa_mapf(jact, iact) = npair
                 rotaa_mapb(npair, 1) = iact ! greater
                 rotaa_mapb(npair, 2) = jact ! lesser
              end if
           end do
        end do
     end do
  end do

  if (iprint > 0) then
     write(6, "(' # ORMAS: active-active rotations: ', i10)") nrotaa
     do irot = 1, nrotaa
        iact = rotaa_mapb(irot, 1)
        jact = rotaa_mapb(irot, 2)
        check = rotaa_mapf(jact, iact)
        if (check /= irot) stop 'inconsistencyr in rotaa_mapf and rotaa_mapb.'
        write(6, "(' irot : ', i5, ' pair: ', 2i5)") irot, iact, jact
     end do
  end if
  
end subroutine ormas_rotaa
!################################################################################
