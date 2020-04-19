!################################################################################
subroutine ormas_rotca()

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : iprint
  use mod_ormas, only : nfcore, ncore, nocc, nrotca, rotca_mapf, rotca_mapb

  implicit none
  integer(c_int) :: ifun, tfun, irot, npair, check

  ! number of core-active rotatins
  nrotca = 0
  do ifun = nfcore + 1, ncore
     do tfun = ncore + 1, nocc
        nrotca = nrotca + 1
     end do
  end do

  allocate(rotca_mapf(1:nocc, 1:nocc))
  allocate(rotca_mapb(1:nrotca, 1:2))

  npair = 0
  do ifun = nfcore + 1, ncore
     do tfun = ncore + 1, nocc
        npair = npair + 1
        rotca_mapf(tfun, ifun) = npair
        rotca_mapb(npair, 1) = tfun
        rotca_mapb(npair, 2) = ifun
     end do
  end do

  if (iprint > 0) then
     write(6, "(' # ORMAS: core-active rotations: ', i10)") nrotca
     do irot = 1, nrotca
        tfun = rotca_mapb(irot, 1)
        ifun = rotca_mapb(irot, 2)
        check = rotca_mapf(tfun, ifun)
        if (check /= irot) stop 'inconsistencyr in rotca_mapf and rotca_mapb.'
        write(6, "(' irot : ', i5, ' pair: ', 2i5)") irot, tfun, ifun
     end do
  end if
  
end subroutine ormas_rotca
!################################################################################
