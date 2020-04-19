!///////////////////////////////////////////////////////////////////////
subroutine wfn_projg(porb, orb)

  use, intrinsic :: iso_c_binding

  use mod_bas, only : ngrid, wgt
  use mod_ormas, only : nocc, nfun

  implicit none
  complex(c_double_complex), intent(in) :: porb(1:ngrid, 1:nfun)
  complex(c_double_complex), intent(inout) :: orb(1:ngrid, 1:nfun)

  integer, external :: omp_get_num_threads
  integer(c_long), external :: util_omp_iproc
  integer(c_long) :: nproc, iproc, llg, ulg

  complex(c_double_complex) :: tmp1, fac
  complex(c_double_complex), allocatable :: s2(:,:,:)
  integer(c_long) :: ifun, jfun, igrid
  integer(c_long), external :: util_omp_nproc

  nproc = util_omp_nproc()
  allocate(s2(1:nfun, 1:nfun, 0:(nproc-1)))

  !$omp parallel default(shared) private(iproc, llg, ulg, tmp1)
  !###########################
  iproc = util_omp_iproc()
  call util_omp_disp(1, ngrid, llg, ulg)

  s2(1:nfun, 1:nfun, iproc) = (0.d+0, 0.d+0)
  do ifun = 1, nfun
     do jfun = 1, nfun
        tmp1 = (0.d+0, 0.d+0)
        do igrid = llg, ulg
           tmp1 = tmp1 + conjg(porb(igrid, jfun)) * orb(igrid, ifun) * wgt(igrid)
        end do
        s2(jfun, ifun, iproc) = tmp1
     end do
  end do
  !###########################
  !$omp end parallel

  do iproc = 1, nproc - 1
     do ifun = 1, nfun
        do jfun = 1, nfun
           s2(jfun, ifun, 0) = s2(jfun, ifun, 0) + s2(jfun, ifun, iproc)
        end do
     end do
  end do

  !$omp parallel default(shared) private(iproc, llg, ulg, fac)
  !###########################
  iproc = util_omp_iproc()
  call util_omp_disp(1, ngrid, llg, ulg)

  do ifun = 1, nocc
     do jfun = 1, nocc
        fac = - s2(jfun, ifun, 0)
        do igrid = llg, ulg
           orb(igrid, ifun) = &
           orb(igrid, ifun) + porb(igrid, jfun) * fac
        end do
     end do
  end do
  !###########################
  !$omp end parallel

  deallocate(s2)

end subroutine wfn_projg
!///////////////////////////////////////////////////////////////////////
