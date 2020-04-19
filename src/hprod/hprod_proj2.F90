!///////////////////////////////////////////////////////////////////////
subroutine hprod_proj2(mpi_gll, mpi_gul, tmat, orb1, orb2)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nfun

  implicit none
  integer(c_int), intent(in) :: mpi_gll, mpi_gul
  complex(c_double_complex), intent(in) :: tmat(1:nfun, 1:nfun)
  complex(c_double_complex), intent(in) :: orb1(mpi_gll:mpi_gul, 1:nfun)
  complex(c_double_complex), intent(inout) :: orb2(mpi_gll:mpi_gul, 1:nfun)

  integer(c_int) :: ifun, jfun, igrid, llg, ulg

  !$omp parallel default(shared) private(llg, ulg)
  !###########################
  call util_omp_disp(mpi_gll, mpi_gul, llg, ulg)

  do ifun = 1, nfun
     do jfun = 1, nfun
        do igrid = llg, ulg
           orb2(igrid, ifun) = orb2(igrid, ifun) &
                             - orb1(igrid, jfun) * tmat(jfun, ifun)
        end do
     end do
  end do
  !###########################
  !$omp end parallel

end subroutine hprod_proj2
!///////////////////////////////////////////////////////////////////////
