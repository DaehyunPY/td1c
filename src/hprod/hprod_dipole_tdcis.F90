!######################################################################
! test FC use nfcore_tdcis instead of nfcore
!######################################################################
subroutine hprod_dipole_tdcis(lfield, dotorb, ci0, dotcic, dip, dip2, vel, vel2)
  !dipole moment in tdcis
  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad, xrad, ecs_flag, irad_ecs
  use mod_sph, only : lmax1
  use mod_ormas, only : nfun
  use mod_hprod, only : orb, orb0rot
  use mod_const, only : runit
  use mod_control, only : tdcis_rvg

  implicit none
  real(c_double), intent(in) :: lfield(1:3, 1:3)
  complex(c_double_complex), intent(in) :: dotorb(1:*)
  complex(c_double_complex), intent(in):: ci0, dotcic
  complex(c_double_complex), intent(out) :: dip, dip2
  complex(c_double_complex), intent(out) :: vel, vel2

! tdcis_sato
!  complex(c_double_complex) :: dmat(1:nfun,1:nfun)
!  ! r|>
!  complex(c_double_complex) :: zphi(1:(nrad - 1), 0:lmax1, 1:nfun)
!  complex(c_double_complex) :: zchi(1:(nrad - 1), 0:lmax1, 1:nfun)
!  complex(c_double_complex) :: dotzchi(1:(nrad - 1), 0:lmax1, 1:nfun)
!  ! p|>
!  complex(c_double_complex) :: pphi(1:(nrad - 1), 0:lmax1, 1:nfun)
!  complex(c_double_complex) :: pchi(1:(nrad - 1), 0:lmax1, 1:nfun)
!  complex(c_double_complex) :: dotpchi(1:(nrad - 1), 0:lmax1, 1:nfun)
  real(c_double) :: radmax
  complex(c_double_complex), allocatable :: dmat(:,:)
  ! r|>
  complex(c_double_complex), allocatable :: zphi(:,:,:)
  complex(c_double_complex), allocatable :: zchi(:,:,:)
  complex(c_double_complex), allocatable :: dotzchi(:,:,:)
  ! p|>
  complex(c_double_complex), allocatable :: pphi(:,:,:)
  complex(c_double_complex), allocatable :: pchi(:,:,:)
  complex(c_double_complex), allocatable :: dotpchi(:,:,:)
  allocate(dmat(1:nfun,1:nfun))
  allocate(zphi(1:(nrad - 1), 0:lmax1, 1:nfun))
  allocate(zchi(1:(nrad - 1), 0:lmax1, 1:nfun))
  allocate(dotzchi(1:(nrad - 1), 0:lmax1, 1:nfun))
  allocate(pphi(1:(nrad - 1), 0:lmax1, 1:nfun))
  allocate(pchi(1:(nrad - 1), 0:lmax1, 1:nfun))
  allocate(dotpchi(1:(nrad - 1), 0:lmax1, 1:nfun))
! tdcis_sato

!write(6,"('for debug F3.0')")
!stop

!tdcis_sato
  if (ecs_flag == 0) then
     radmax = xrad(nrad)
  else 
     radmax = xrad(irad_ecs-1)
  end if
!tdcis_sato

  ! complex(c_double_complex) :: pos(1:nfun,1:nfun), mom(1:nfun,1:nfun)
  call hprod_mkovlp(radmax, orb, dotorb, dmat)

  !  call hprod_popgs_tdcis(orb)

!write(6,"('for debug F3.1')")
!stop
  ! z|>
  dip = 0d0
  dip2 = 0d0
  zphi = 0d0
  zchi = 0d0
  dotzchi = 0d0
  call hprod_zprod_all(runit, orb0rot, zphi)
  call hprod_zprod_all(runit, orb, zchi)

  ! p|>
  vel = 0d0
  vel2 = 0d0
  pphi = 0d0
  pchi = 0d0
  dotpchi = 0d0
  call hprod_pzprod_all(runit, orb0rot, pphi)
  call hprod_pzprod_all(runit, orb, pchi)

!write(6,"('for debug F3.2')")
!stop
  ! <z>
  ! write(6, "('hprod_tdcis_r')")
  call hprod_op1e_tdcis(orb0rot, orb, zphi, zchi, ci0, dip)
  call hprod_dotop1e_tdcis(orb0rot, orb, zphi, zchi, dotorb, ci0, dotcic, dmat, dip2) 

!write(6,"('for debug F3.3')")
!stop
  ! <p>
  ! write(6, "('hprod_tdcis_p')")
  call hprod_op1e_tdcis(orb0rot, orb, pphi, pchi, ci0, vel) ! kinetic momentum, not canonical momentum
  call hprod_dotop1e_tdcis(orb0rot, orb, pphi, pchi, dotorb, ci0, dotcic, dmat, vel2) 

!write(6,"('for debug F3.3')")
!stop
  ! correction for rVG
  if(tdcis_rvg) then
     call hprod_vel_add_tdcis(lfield, zphi, zchi, ci0, dip2)
     call hprod_acc_add_tdcis(lfield, zphi, pchi, ci0, vel2)
  end if

! tdcis_sato
  deallocate(dmat)
  deallocate(zphi)
  deallocate(zchi)
  deallocate(dotzchi)
  deallocate(pphi)
  deallocate(pchi)
  deallocate(dotpchi)
! tdcis_sato

!write(6,"('for debug F3.4')")
!stop

end subroutine hprod_dipole_tdcis
!######################################################################
complex(c_double_complex) function hprod_trace_tdcis(mat)

  use, intrinsic :: iso_c_binding
  use mod_const, only : czero
  use mod_ormas, only : nfun, nfcore 
  implicit none
  complex(c_double_complex), intent(in) :: mat(1:nfun, 1:nfun)

  integer(c_int) :: ifun
  complex(c_double_complex) :: tval

  tval = czero
  ! do ifun = 1, nfun
  do ifun = nfcore + 1, nfun
     tval = tval + mat(ifun, ifun)
  end do

  hprod_trace_tdcis = tval
  return

end function hprod_trace_tdcis
!######################################################################
complex(c_double_complex) function hprod_trace2_tdcis(rmax, wfn, hwfn)

  use, intrinsic :: iso_c_binding
  use mod_const, only : czero
  use mod_ormas, only : nfun, nfcore 
  implicit none
  real(c_double), intent(in) :: rmax
  complex(c_double_complex), intent(in) :: wfn(1:*), hwfn(1:*)

  integer(c_int) :: ifun
  complex(c_double_complex) :: mat(1:nfun, 1:nfun)
  complex(c_double_complex) :: tval

  tval = czero
  call hprod_mkovlp(rmax, wfn, hwfn, mat)

  ! do ifun = 1, nfun
  do ifun = nfcore + 1, nfun
     tval = tval + mat(ifun, ifun)
  end do

  hprod_trace2_tdcis = tval
  return

end function hprod_trace2_tdcis
!######################################################################
