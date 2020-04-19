!######################################################################
real(c_double) function hprod_ene_fcore(doall, orb, h0orb, h1orb, gorb)

  use, intrinsic :: iso_c_binding
  use mod_const, only : zero, czero
  use mod_bas, only : nbas
  use mod_control, only : ioorot
  use mod_ormas, only : nfcore2, nfcore, nfun
  use mod_rad, only : ecs_flag

  implicit none
  logical(c_bool), intent(in) :: doall
  complex(c_double_complex), intent(in) ::   orb(1:nbas, 1:nfun)
  complex(c_double_complex), intent(in) :: h0orb(1:nbas, 1:nfun)
  complex(c_double_complex), intent(in) :: h1orb(1:nbas, 1:nfun)
  complex(c_double_complex), intent(in) ::  gorb(1:nbas, 1:nfun)

  real(c_double), external :: hprod_ene_fcx
  integer(c_int), external :: util_omp_nproc
  integer(c_int), external :: util_omp_iproc
  complex(c_double_complex) :: tmp, h01, h012, enep
  integer(c_int) :: nproc, iproc, ifun, ibas, llb, ulb

!nyi  real(c_double), external :: hprod_ene_fcore_ecs
!nyi  if (ecs_flag == 1) then
!nyi     hprod_ene_fcore = hprod_ene_fcore_ecs(orb, h0orb, h1orb, gorb)
!nyi     return
!nyi  end if

!debug
!  do ifun = 1, nfun
!     write(6, "('# ifun = ', i5)") ifun
!     do ibas = 1, nbas
!        write(6, "(i5, 6f20.10)") ibas, h0orb(ibas, ifun), h1orb(ibas, ifun), gorb(ibas, ifun)
!     end do
!  end do
!  stop
!debug

  if (nfcore == 0) then
     hprod_ene_fcore = zero
     return
  end if

  nproc = util_omp_nproc()
  enep = czero

  !$omp parallel default(shared) private(iproc, llb, ulb, tmp, h01, h012) reduction(+:enep)
  !###########################
  iproc = util_omp_iproc()
  call util_omp_disp(1, nbas, llb, ulb)
  tmp = czero
  do ifun = 1, nfcore
     if (doall .or. ioorot == 0) then
        ! <j|h0+h1+h2|i>
        do ibas = llb, ulb
           h01 = h1orb(ibas, ifun) + h0orb(ibas, ifun)
           h012 = gorb(ibas, ifun) + h01 + h01
           tmp = tmp + conjg(orb(ibas, ifun)) * h012
        end do
     else if (ioorot == 2) then
        ! <j|h1+h2|i>
        do ibas = llb, ulb
           h01 = h1orb(ibas, ifun)
           h012 = gorb(ibas, ifun) + h01 + h01
           tmp = tmp + conjg(orb(ibas, ifun)) * h012
        end do
     else if (ioorot == 1) then
        ! <j|h2|i>
        do ibas = llb, ulb
           h012 = gorb(ibas, ifun)
           tmp = tmp + conjg(orb(ibas, ifun)) * h012
        end do
     end if
  end do
  enep = enep + tmp
  !###########################
  !$omp end parallel

!DEBUG
! write(6, "('ene_fcore_1 = ', f20.10)") dble(enep)
! write(6, "('hprod_ene_fcore: skip hprod_ene_fcx')") dble(enep)
!  if (nfcore2 > 0) then
!     enep = enep + hprod_ene_fcx()
!  end if
!DEBUG

  hprod_ene_fcore = dble(enep)
  return

end function hprod_ene_fcore
!######################################################################
