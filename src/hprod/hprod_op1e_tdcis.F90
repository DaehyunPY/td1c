!######################################################################
! test FC use nfcore_tdcis instead of nfcore
!######################################################################
subroutine hprod_op1e_tdcis_old(phi, chi, ophi, ochi, ci0, dip)
  !<o> in tdcis general
  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad, xrad, ecs_flag, irad_ecs
  use mod_sph, only : lmax1
  use mod_ormas, only : nfun, ncore, nact
  use mod_hprod, only : ridm_tdcis
  use mod_const, only : czero

  implicit none
  complex(c_double_complex), intent(in) :: phi(1:(nrad - 1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: chi(1:(nrad - 1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: ophi(1:(nrad - 1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: ochi(1:(nrad - 1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: ci0
  complex(c_double_complex), external :: hprod_trace_tdcis
  complex(c_double_complex), intent(inout) :: dip

  integer(c_int) :: ifun, jfun
  complex(c_double_complex) :: ophiphi(1:nfun, 1:nfun)
  complex(c_double_complex) :: ochiphi(1:nfun, 1:nfun)
  complex(c_double_complex) :: ochichi(1:nfun, 1:nfun)
  real(c_double) :: radmax
  complex(c_double_complex) :: sum1, sum2, sum3, sum4, sum5
  complex(c_double_complex) :: tmp1, tmp2, tmp3

  ! initialization
  sum1 = czero
  sum2 = czero
  sum3 = czero
  sum4 = czero
  sum5 = czero
  tmp1 = 2d0*abs(ci0)*abs(ci0)
  tmp2 = 2d0*sqrt(2d0)
  tmp3 = czero
  dip = czero

!tdcis_sato
  if (ecs_flag == 0) then
     radmax = xrad(nrad)
  else 
     radmax = xrad(irad_ecs-1)
  end if
!tdcis_sato

  ! matrix
  call hprod_mkovlp(radmax, phi, ophi, ophiphi)
  call hprod_mkovlp(radmax, chi, ochi, ochichi)
  call hprod_mkovlp(radmax, chi, ophi, ochiphi)

  dip = hprod_trace_tdcis(ophiphi)*2d0*abs(ci0)*abs(ci0) + hprod_trace_tdcis(ochichi) &
       + real(hprod_trace_tdcis(ochiphi)*ci0) + hprod_trace_tdcis(ridm_tdcis)*hprod_trace_tdcis(ophiphi)
  do ifun = 1, nfun
     do jfun = 1, nfun
        sum4 = sum4 + ridm_tdcis(jfun, ifun)*ophiphi(ifun,jfun)
     end do
  end do
  dip = dip - sum4
  ! write(6 , "('hprod_dip_tdcis:debug', 6es20.10)") ci0, ridm_tdcis(1,1), dip

end subroutine hprod_op1e_tdcis_old
!######################################################################
subroutine hprod_op1e_tdcis(phi, chi, ophi, ochi, ci0, dip)
  !<o> in tdcis <phi|o|phi> = 0
  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad, xrad, ecs_flag, irad_ecs
  use mod_sph, only : lmax1
  use mod_ormas, only : nfun, nfcore
  use mod_hprod, only : ridm_tdcis
  use mod_const, only : czero

  implicit none
  complex(c_double_complex), intent(in) :: phi(1:(nrad - 1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: chi(1:(nrad - 1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: ophi(1:(nrad - 1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: ochi(1:(nrad - 1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: ci0
  complex(c_double_complex), external :: hprod_trace2_tdcis
  complex(c_double_complex), intent(out) :: dip

  integer(c_int) :: ifun, jfun
  complex(c_double_complex) :: sum
  complex(c_double_complex) :: ophiphi(1:nfun, 1:nfun)
  complex(c_double_complex) :: sum1, sum2
!tdcis_sato
  real(c_double) :: radmax
  if (ecs_flag == 0) then
     radmax = xrad(nrad)
  else 
     radmax = xrad(irad_ecs-1)
  end if
!tdcis_sato

  sum1 = czero
  sum2 = czero
  ophiphi = czero
  sum = czero

  dip = hprod_trace2_tdcis(radmax, chi, ochi) + 2d0*sqrt(2d0)*real(ci0*hprod_trace2_tdcis(radmax, chi, ophi))

  call hprod_mkovlp(radmax, phi, ophi, ophiphi)
  do ifun = 1, nfun
     do jfun = 1, nfun
        sum = sum + ridm_tdcis(ifun, jfun)*ophiphi(jfun, ifun)
     end do
  end do

  dip = dip - sum

  !test for VG
  sum1 = hprod_trace2_tdcis(radmax, chi, ochi) 
  sum2 = 2d0*sqrt(2d0)*real(ci0*hprod_trace2_tdcis(radmax, chi, ophi))
  ! write(6, "('hprod_op1e_tdcis:debug', 4e20.10)") sum1, sum2
  !test for VG

end subroutine hprod_op1e_tdcis
!######################################################################
subroutine hprod_dotop1e_tdcis(phi, chi, ophi, ochi, dotchi, ci0, dc, dmat, vel)
  !<o> in tdcis
  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad, xrad, ecs_flag, irad_ecs
  use mod_sph, only : lmax1
  use mod_ormas, only : nfun, nfcore
  use mod_hprod, only : ridm_tdcis
  use mod_const, only : czero

  implicit none
  complex(c_double_complex), intent(in) :: phi(1:(nrad - 1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: chi(1:(nrad - 1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: ophi(1:(nrad - 1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: ochi(1:(nrad - 1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: dotchi(1:(nrad - 1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: ci0, dc
  complex(c_double_complex), intent(in) :: dmat(1:nfun, 1:nfun)
  complex(c_double_complex), external :: hprod_trace2_tdcis
  complex(c_double_complex), intent(out) :: vel

  integer(c_int) :: ifun, jfun
  complex(c_double_complex) :: ophiphi(1:nfun, 1:nfun)
  complex(c_double_complex) :: sum1, sum2

  complex(c_double_complex) :: t1, t2, t3
!tdcis_sato
  real(c_double) :: radmax
  if (ecs_flag == 0) then
     radmax = xrad(nrad)
  else 
     radmax = xrad(irad_ecs-1)
  end if
!tdcis_sato

  ! initialization
  sum1 = czero
  sum2 = czero
  vel = czero
  t1 = czero
  t2 = czero
  t3 = czero

  ! test
  t1 = 2d0*hprod_trace2_tdcis(radmax, ochi, dotchi)
  t2 = 2d0*sqrt(2d0)*dc*hprod_trace2_tdcis(radmax, chi, ophi)
  t3 = 2d0*sqrt(2d0) *ci0* hprod_trace2_tdcis(radmax, dotchi, ophi)

  ! test
  sum1 = hprod_trace2_tdcis(radmax, ochi, dotchi) &
       + sqrt(2d0)*dc*hprod_trace2_tdcis(radmax, chi, ophi) &
       + sqrt(2d0) *ci0* hprod_trace2_tdcis(radmax, dotchi, ophi)

  ! matrix
  call hprod_mkovlp(radmax, phi, ophi, ophiphi)
  ! do jfun = 1, nfun
  !    do ifun = 1, nfun
  do jfun = nfcore + 1, nfun
     do ifun = nfcore + 1, nfun
        sum2 = sum2 + dmat(jfun, ifun)*ophiphi(ifun,jfun)
     end do
  end do
  ! write(6, "('hprod_dotop1e_tdcis:debug', 4e20.10)") real(t1), real(t2), real(t3), real(sum2)
  vel = 2d0*real(sum1 - sum2)

end subroutine hprod_dotop1e_tdcis
!######################################################################
subroutine hprod_vel_add_tdcis(lfield, zphi, zchi, ci0, vel)
  ! <v> correcton
  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad, xrad, ecs_flag, irad_ecs
  use mod_sph, only : lmax1
  use mod_ormas, only : nfun, nfcore
  use mod_hprod, only : ridm_tdcis
  use mod_const, only : czero

  implicit none
  real(c_double), intent(in) :: lfield(1:3, 1:3)
  complex(c_double_complex), intent(in) :: zphi(1:(nrad - 1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: zchi(1:(nrad - 1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: ci0
  complex(c_double_complex), external :: hprod_trace2_tdcis
  complex(c_double_complex), intent(inout) :: vel

  integer(c_int) :: ifun
  complex(c_double_complex) :: cval
  complex(c_double_complex) :: efield, afield
  complex(c_double_complex) :: mat1(1:nfun, 1:nfun)
  complex(c_double_complex) :: mat2(1:nfun, 1:nfun)
!tdcis_sato
  real(c_double) :: radmax
  if (ecs_flag == 0) then
     radmax = xrad(nrad)
  else 
     radmax = xrad(irad_ecs-1)
  end if
!tdcis_sato

  efield = lfield(3, 2)
  cval = hprod_trace2_tdcis(radmax, zchi, zphi)

  ! cval = czero
  ! do ifun = 1, nfun 
  !    cval = cval + mat1(ifun, ifun) 
  ! end do
  cval = -2d0*sqrt(2d0)*efield*aimag(ci0*cval)
  vel = vel + cval

end subroutine hprod_vel_add_tdcis
!######################################################################
subroutine hprod_acc_add_tdcis(lfield, zphi, pchi, ci0, acc)
  ! <a> correcton
  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad, xrad, ecs_flag, irad_ecs
  use mod_sph, only : lmax1
  use mod_ormas, only : nfun, nfcore
  use mod_hprod, only : ridm_tdcis
  use mod_const, only : czero

  implicit none
  real(c_double), intent(in) :: lfield(1:3, 1:3)
  complex(c_double_complex), intent(in) :: zphi(1:(nrad - 1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: pchi(1:(nrad - 1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: ci0
  complex(c_double_complex), external :: hprod_trace_tdcis
  complex(c_double_complex), external :: hprod_trace2_tdcis
  complex(c_double_complex), intent(inout) :: acc

  integer(c_int) :: ifun
  complex(c_double_complex) :: cval, tval
  complex(c_double_complex) :: efield, afield
  complex(c_double_complex) :: mat1(1:nfun, 1:nfun)
!tdcis_sato
  real(c_double) :: radmax
  if (ecs_flag == 0) then
     radmax = xrad(nrad)
  else 
     radmax = xrad(irad_ecs-1)
  end if
!tdcis_sato

  efield = lfield(3, 2)
  
  cval = hprod_trace2_tdcis(radmax, pchi, zphi)
  tval = hprod_trace_tdcis(ridm_tdcis)

  ! cval = czero
  ! tval = czero
  ! call hprod_mkovlp(radmax, pchi, zphi, mat1) !<chi|p z|phi>

  ! do ifun = 1, nfun
  !    cval = cval + mat1(ifun, ifun)
  !    tval = tval + ridm_tdcis(ifun, ifun)
  ! end do
  cval = -2d0*sqrt(2d0)*efield*aimag(ci0*cval)
  cval = cval - tval * efield 
  acc = acc + cval

end subroutine hprod_acc_add_tdcis
!######################################################################
