!///////////////////////////////////////////////////////////////////////
subroutine hprod_proj(tmat, orb1, orb2)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : nbas, mval
  use mod_ormas, only : nfun

  implicit none
  ! ### args ###
  complex(c_double_complex), intent(in) :: tmat(1:nfun, 1:nfun)
  complex(c_double_complex), intent(in) :: orb1(1:nbas, 1:nfun)
  complex(c_double_complex), intent(inout) :: orb2(1:nbas, 1:nfun)
  ! ### local ###
  integer(c_long) :: ifun, jfun, ibas, mi, mj, m, llb, ulb

  !$omp parallel default(shared) private(mi, mj, m, llb, ulb)
  !###########################
  call util_omp_disp(1, nbas, llb, ulb)
  do ifun = 1, nfun
     mi = mval(ifun)

     do jfun = 1, nfun
        mj = mval(jfun)

        if (mj == mi) then
           do ibas = llb, ulb
              orb2(ibas, ifun) = orb2(ibas, ifun) &
                               - orb1(ibas, jfun) * tmat(jfun, ifun)
           end do
        end if
     end do
  end do
  !###########################
  !$omp end parallel

end subroutine hprod_proj
!///////////////////////////////////////////////////////////////////////
