!######################################################################
subroutine hprod_mkovlp2(mpi_gll, mpi_gul, wfn, hwfn, ovlp)

  use, intrinsic :: iso_c_binding
  use mod_const, only : czero
  use mod_ormas, only : nfun

  implicit none
  integer(c_int), intent(in) :: mpi_gll, mpi_gul
  complex(c_double_complex), intent(in) :: wfn(mpi_gll:mpi_gul, 1:nfun)
  complex(c_double_complex), intent(in) :: hwfn(mpi_gll:mpi_gul, 1:nfun)
  complex(c_double_complex), intent(out) :: ovlp(1:nfun, 1:nfun)

  integer(c_int), external :: util_omp_iproc
  integer(c_int), external :: util_omp_nproc
  integer(c_int) :: ifun, jfun, iproc, nproc, ngrid, llg, ulg
  complex(c_double_complex), allocatable :: ovlpp(:,:,:)

  nproc = util_omp_nproc()
  ngrid = mpi_gul - mpi_gll + 1

  allocate(ovlpp(1:nfun, 1:nfun, 0:(nproc-1)))

  !$omp parallel default(shared) private(iproc, llg, ulg)
  !###########################
  iproc = util_omp_iproc()
  call util_omp_disp(mpi_gll, mpi_gul, llg, ulg)

  ovlpp(1:nfun, 1:nfun, iproc) = czero
  call hprod_mkovlp2p(mpi_gll, mpi_gul, wfn, hwfn, ovlpp(1,1,iproc), llg, ulg)
  !###########################
  !$omp end parallel

  ovlp(1:nfun, 1:nfun) = czero
  do iproc = 0, nproc - 1
     do ifun = 1, nfun
        do jfun = 1, nfun
           ovlp (jfun, ifun) = &
         & ovlp (jfun, ifun) + &
         & ovlpp(jfun, ifun, iproc)
        end do
     end do
  end do

  deallocate(ovlpp)

end subroutine hprod_mkovlp2
!######################################################################
subroutine hprod_mkovlp2p(mpi_gll, mpi_gul, wfn, hwfn, ovlp, llg, ulg)

  use, intrinsic :: iso_c_binding
  use mod_const, only : czero
  use mod_ormas, only : nfun
  use mod_bas, only : wgt

  implicit none
  integer(c_int), intent(in) :: mpi_gll, mpi_gul, llg, ulg
  complex(c_double_complex), intent(in) :: wfn(mpi_gll:mpi_gul, 1:nfun)
  complex(c_double_complex), intent(in) :: hwfn(mpi_gll:mpi_gul, 1:nfun)
  complex(c_double_complex), intent(inout) :: ovlp(1:nfun, 1:nfun)

  integer(c_int) :: ifun, jfun, igrid

  do ifun = 1, nfun
     do jfun = 1, nfun
        do igrid = llg, ulg
           ovlp(jfun, ifun) = ovlp(jfun, ifun) + conjg(wfn(igrid, jfun)) * hwfn(igrid, ifun) * wgt(igrid)
        end do
     end do
  end do

end subroutine hprod_mkovlp2p
!######################################################################
