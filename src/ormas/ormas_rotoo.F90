!################################################################################
subroutine ormas_rotoo()

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : iprint
  use mod_ormas, only : ncore, nfun
  use mod_ormas, only : nrotoo, rotoo_mapf, rotoo_mapb
  use mod_ormas, only : nrotca, rotca_mapb
  use mod_ormas, only : nrotaa, rotaa_mapb

  implicit none
  integer(c_int) :: irot, ifun, jfun, check

  call ormas_rotca()
  call ormas_rotaa()

  nrotoo = nrotca + nrotaa
  allocate(rotoo_mapf(1:nfun, 1:nfun))
  allocate(rotoo_mapb(1:nrotoo, 1:2))

  do irot = 1, nrotca
     ifun = rotca_mapb(irot, 1)
     jfun = rotca_mapb(irot, 2)
     rotoo_mapf(jfun, ifun) = irot
     rotoo_mapb(irot, 1) = ifun
     rotoo_mapb(irot, 2) = jfun
  end do

  do irot = 1, nrotaa
     ifun = ncore + rotaa_mapb(irot, 1)
     jfun = ncore + rotaa_mapb(irot, 2)
     rotoo_mapf(jfun, ifun) = nrotca + irot
     rotoo_mapb(nrotca + irot, 1) = ifun
     rotoo_mapb(nrotca + irot, 2) = jfun
  end do

  if (iprint > 1) then
     write(6, "('# ORMAS: occ-occ rotations: ', i10)") nrotoo
     do irot = 1, nrotoo
        ifun = rotoo_mapb(irot, 1)
        jfun = rotoo_mapb(irot, 2)
        check = rotoo_mapf(jfun, ifun)
        if (check /= irot) stop 'inconsistencyr in rotoo_mapf and rotoo_mapb.'
        write(6, "(' irot : ', i5, ' pair: ', 2i5)") irot, ifun, jfun
     end do
  end if

end subroutine ormas_rotoo
!################################################################################
