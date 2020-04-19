!######################################################################
subroutine bas_pp_genkb(orb)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad, xrad, wrad
  use mod_sph, only : lmax1
  use mod_bas, only : lval, mval, pp_vlocHF
  use mod_ormas, only : nfun
  use mod_const, only : czero

  implicit none
  !--------------------------------------------------------------------
  complex(c_double_complex), intent(in) :: orb(1:(nrad-1), 0:lmax1, 1:nfun)
  !--------------------------------------------------------------------
  integer(c_long) :: ifun, irad, myl,mym
  complex(c_double_complex) :: norm(1:nfun), ovlp(1:nfun)
  complex(c_double_complex), allocatable :: chi(:,:)

  allocate(chi(1:(nrad-1), 1:nfun))

  chi = 0d0
  norm = 0d0
  ovlp = 0d0

  do ifun = 1, nfun
  do irad = 1, nrad-1
     myl = lval(ifun)
     mym = mval(ifun)
     chi(irad,ifun) = orb(irad,myl,ifun)*pp_vlocHF(irad,min(myl,2))
     norm(ifun) = norm(ifun) + conjg(orb(irad,myl,ifun))*orb(irad,myl,ifun)
     ovlp(ifun) = ovlp(ifun) + conjg(orb(irad,myl,ifun))*chi(irad,ifun)
  end do
  end do

  write(6,"('ProjKB-norm: ')",advance='no')
  do ifun = 1, nfun
     write(6,"(f20.10)",advance='no') dble(norm(ifun))
  end do
  write(6,*)
  write(6,"('ProjKB-ovlp: ')",advance='no')
  do ifun = 1, nfun
     write(6,"(f20.10)",advance='no') dble(ovlp(ifun))
  end do
  write(6,*)

  do irad = 1, nrad-1
     write(6,"('ProjKB-chi: ',f20.10)",advance='no') xrad(irad)
     do ifun = 1, nfun
        write(6,"(2f20.10)",advance='no') chi(irad,ifun)/sqrt(wrad(irad))
     end do
     write(6,*)
  end do

  deallocate(chi)

end subroutine bas_pp_genkb
!######################################################################
