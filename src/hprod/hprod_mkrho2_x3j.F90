!######################################################################
subroutine hprod_mkrho2_x3j()

  use, intrinsic :: iso_c_binding
  use mod_hprod, only : orb, rho2

  implicit none
  call hprod_mkrho2p_x3j('tot', orb, rho2)

end subroutine hprod_mkrho2_x3j
!######################################################################
subroutine hprod_mkrho2_x3j_dyn()

  use, intrinsic :: iso_c_binding
  use mod_hprod, only : orb, rho2

  implicit none

  call hprod_mkrho2p_x3j('dyn', orb, rho2)
  call hprod_mkrho2p_x3j('fc1dyn', orb, rho2)

!debug
!  write(6, "('hprod_mkrho2_x3j_dyn: rho2')")
!  call hprod_printrho(rho2)
!  stop
!debug

end subroutine hprod_mkrho2_x3j_dyn
!######################################################################
subroutine hprod_mkrho2p_x3j(v2_type, orb, rho2)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : mval
  use mod_sph, only : lmax1,lmax2,mmax1,mmax2,sph_gaunt
  use mod_const, only : ctwo, pi
  use mod_ormas, only : nfun, nfcore, nfcore1, nfcore2
  use mod_rad, only : nrad, nradfc

  implicit none
  character(len=*), intent(in) :: v2_type
  complex(c_double_complex), intent(in) :: orb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(out) :: rho2(1:(nrad-1), 0:lmax2, 1:nfun, 1:nfun)
  integer(c_long) :: dim,ifun_ll,ifun_ul,jfun_ll,jfun_ul,llr,ulr,ifun,jfun,irad
  integer(c_long) :: li,lj,lij,mij,mfac,lllij,ullij
  complex(c_double_complex) :: rhoji(1:(nrad-1))

  if (trim(v2_type) == 'tot') then
     dim = nrad - 1
     ifun_ll = 1
     ifun_ul = nfun
     jfun_ll = 1
     jfun_ul = nfun
  else if (trim(v2_type) == 'fcfc') then
     dim = nradfc
!    dim = nrad - 1
     ifun_ll = 1
     ifun_ul = nfcore
     jfun_ll = 1
     jfun_ul = nfcore
  else if (trim(v2_type) == 'fc1fc1') then
     dim = nradfc
!    dim = nrad - 1
     ifun_ll = nfcore2 + 1
     ifun_ul = nfcore
     jfun_ll = nfcore2 + 1
     jfun_ul = nfcore
  else if (trim(v2_type) == 'fc2fc2') then
     dim = nradfc
!    dim = nrad - 1
     ifun_ll = 1
     ifun_ul = nfcore2
     jfun_ll = 1
     jfun_ul = nfcore2
  else if (trim(v2_type) == 'fc1dyn') then
     dim = nradfc
     ifun_ll = nfcore2 + 1
     ifun_ul = nfcore
     jfun_ll = nfcore + 1
     jfun_ul = nfun
  else if (trim(v2_type) == 'fc2dyn') then
     dim = nradfc
     ifun_ll = 1
     ifun_ul = nfcore2
     jfun_ll = nfcore + 1
     jfun_ul = nfun
  else if (trim(v2_type) == 'dyn') then
     dim = nrad - 1
     ifun_ll = nfcore + 1
     ifun_ul = nfun
     jfun_ll = nfcore + 1
     jfun_ul = nfun
  else
     stop "bad v2_type in hprod_mkrho2p"
  end if

  !$omp parallel default(shared) private(llr,ulr,mij,mfac,lllij,ullij)
  !###########################
  call util_omp_disp(1, dim, llr, ulr)
  do ifun = ifun_ll, ifun_ul
     do jfun = max(jfun_ll, ifun), max(jfun_ul, ifun)
        mij = mval(ifun)-mval(jfun)
        mfac = (-1)**mij
        if (abs(mij) > mmax2) cycle
        do li = abs(mval(ifun)), lmax1
           do lj = abs(mval(jfun)), lmax1
              lllij = max(abs(mij),abs(li-lj))
              ullij = min(lmax2,li+lj)
              do irad = llr, ulr
                 rhoji(irad) = conjg(orb(irad,lj,jfun))*orb(irad,li,ifun)*mfac
!                 rhoji(irad) = conjg(orb(irad,lj,jfun))*orb(irad,li,ifun)
!                 rhoji(irad) = conjg(orb(irad,lj,jfun))*orb(irad,li,ifun)*sqrt(2*pi)
              end do
              do lij = lllij, ullij
                 if (mod(li+lj+lij,2).ne.0) cycle
                 do irad = llr, ulr
                    rho2(irad,lij,jfun,ifun) = &
                    rho2(irad,lij,jfun,ifun) + rhoji(irad) &
                    * sph_gaunt(lij,lj,li,mval(jfun),mval(ifun))
                 end do
              end do
           end do
        end do
     end do
  end do
  !###########################
  !$omp end parallel

!old  !$omp parallel default(shared) private(llr,ulr,mij,lllij,ullij)
!old  !###########################
!old  call util_omp_disp(1, dim, llr, ulr)
!old  do ifun = ifun_ll, ifun_ul
!old     do jfun = max(jfun_ll, ifun), max(jfun_ul, ifun)
!old!     do jfun = jfun_ll, jfun_ul
!old        mij = mval(ifun)-mval(jfun)
!old        if (abs(mij) > mmax2) cycle
!old        do li = abs(mval(ifun)), lmax1
!old           do lj = abs(mval(jfun)), lmax1
!old              lllij = max(abs(mij),abs(li-lj))
!old              ullij = min(lmax2,li+lj)
!old              do irad = llr, ulr
!old!NOTE HERE. SEE ALSO hprod_mfprodx3j and bas_gen_d2fac
!old                 rhoji(irad) = conjg(orb(irad,lj,jfun))*orb(irad,li,ifun)
!old!                 rhoji(irad) = conjg(orb(irad,lj,jfun))*orb(irad,li,ifun)*sqrt(2*pi)
!old!NOTE HERE. SEE ALSO hprod_mfprodx3j and bas_gen_d2fac
!old              end do
!old              do lij = lllij, ullij
!old                 if (mod(li+lj+lij,2).ne.0) cycle
!old                 do irad = llr, ulr
!old                    rho2(irad,lij,jfun,ifun) = &
!old                    rho2(irad,lij,jfun,ifun) + rhoji(irad) &
!old                    * sph_gaunt(lij,li,lj,mval(ifun),mval(jfun))
!old!the same                    rho2(irad,lij,jfun,ifun) = &
!old!the same                    rho2(irad,lij,jfun,ifun) + rhoji(irad) &
!old!the same                    * sph_gaunt(lij,li,lj,-mval(ifun),-mval(jfun))
!old                 end do
!old              end do
!old           end do
!old        end do
!old     end do
!old  end do
!old  !###########################
!old  !$omp end parallel

end subroutine hprod_mkrho2p_x3j
!######################################################################
