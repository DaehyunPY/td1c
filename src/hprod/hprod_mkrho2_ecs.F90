subroutine hprod_mkrho2_dyn_ecs()

  use, intrinsic :: iso_c_binding
  use mod_hprod, only : orb, orbg, rho2, v2ang
  use mod_hprod, only : orbe, orbo, v2ange, v2ango, rho23j

  implicit none
  call hprod_mkrho2p_ecs('dyn', orbg, v2ang)
  call bas_ang2sph2_dyn_ecs(v2ang, rho2);
  call hprod_mkrho2p_ecs('fc1dyn', orbg, v2ang)
  call bas_ang2sph2_fc1dyn_ecs(v2ang, rho2);

end subroutine hprod_mkrho2_dyn_ecs
!######################################################################
subroutine hprod_mkrho2p_ecs(v2_type, orbg, rho2g)

  use, intrinsic :: iso_c_binding
  use mod_sph, only : nlat
  use mod_const, only : ctwo
  use mod_ormas, only : nfun, nfcore, nfcore1, nfcore2
  use mod_rad, only : xrad, wrad, nrad, nradfc, irad_ecs

  implicit none
  character(len=*), intent(in) :: v2_type
  complex(c_double_complex), intent(in) :: orbg(1:(nrad-1), 1:nlat, 1:nfun)
  complex(c_double_complex), intent(out) :: rho2g(1:(nrad-1), 1:nlat, 1:nfun, 1:nfun)
  integer(c_long) :: dim, ifun_ll, ifun_ul, jfun_ll, jfun_ul, llr, ulr, ifun, jfun, ilat, irad

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

  !$omp parallel default(shared) private(llr, ulr)
  !###########################
  call util_omp_disp(1, irad_ecs-1, llr, ulr)
  do ifun = ifun_ll, ifun_ul
     do jfun = max(jfun_ll, ifun), max(jfun_ul, ifun)
        do ilat = 1, nlat
           do irad = llr, ulr
             rho2g(irad, ilat, jfun, ifun) = conjg(orbg(irad, ilat, jfun)) &
                                                 * orbg(irad, ilat, ifun)
             rho2g(irad, ilat, ifun, jfun) = conjg(orbg(irad, ilat, ifun)) &
                                                 * orbg(irad, ilat, jfun)
           end do
        end do
     end do
  end do
  !###########################
  !$omp end parallel

end subroutine hprod_mkrho2p_ecs
