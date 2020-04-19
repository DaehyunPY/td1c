!///////////////////////////////////////////////////////////////////////
subroutine hprod_mklzmat(orb, int1e)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nfun
  use mod_rad, only : nrad
  use mod_sph, only : lmax1
  use mod_bas, only : mval
  use mod_const, only : half, czero, iunit

  implicit none
  complex(c_double_complex), intent(in) :: orb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(out) :: int1e(1:nfun, 1:nfun)

  integer(c_int), external :: util_omp_iproc, util_omp_nproc
  integer(c_int) :: ifun, jfun, irad, l, m, mi, mj, nproc, iproc, llr, ulr
  complex(c_double_complex) :: tmp1z
  complex(c_double_complex), allocatable :: oazp(:,:,:)

  nproc = util_omp_nproc()
  allocate(oazp(1:nfun, 1:nfun, 0:(nproc-1)))
  oazp(1:nfun, 1:nfun, 0:(nproc-1)) = czero

  !$omp parallel default(shared) private(iproc, llr, ulr, mi, mj, m, tmp1z)
  !###########################
  iproc = util_omp_iproc()
  call util_omp_disp(1, nrad - 1, llr, ulr)

  oazp(1:nfun, 1:nfun, iproc) = czero
  do ifun = 1, nfun
     mi = mval(ifun)
     do jfun = 1, nfun
        mj = mval(jfun)
        if (mi == mj) then
           m = mi
           tmp1z = czero
           do l = abs(m), lmax1
              do irad = llr, ulr
                 tmp1z = tmp1z + conjg(orb(irad, l, jfun)) * orb(irad, l, ifun)
              end do
           end do
           oazp(jfun, ifun, iproc) = dble(m) * tmp1z
        end if
     end do
  end do
  !###########################
  !$omp end parallel

  do iproc = 1, nproc - 1
     do ifun = 1, nfun
        do jfun = 1, nfun
           oazp(jfun, ifun, 0) = oazp(jfun, ifun, 0) + oazp(jfun, ifun, iproc)
        end do
     end do
  end do

  int1e(1:nfun, 1:nfun) = oazp(1:nfun, 1:nfun, 0)

  !DEBUG
!  write(6, "('hprod_mklzmat. oax (R)')")
!  do ifun = 1, nfun
!     do jfun = 1, nfun
!        write(6, "(f20.10)", advance = 'no') dble(oaz(jfun, ifun, 1))
!     end do
!     write(6, *)
!  end do
!  write(6, "('hprod_mklzmat. oax (I)')")
!  do ifun = 1, nfun
!     do jfun = 1, nfun
!        write(6, "(f20.10)", advance = 'no') aimag(oaz(jfun, ifun, 1))
!     end do
!     write(6, *)
!  end do
!  write(6, "('hprod_mklzmat. oay (R)')")
!  do ifun = 1, nfun
!     do jfun = 1, nfun
!        write(6, "(f20.10)", advance = 'no') dble(oaz(jfun, ifun, 2))
!     end do
!     write(6, *)
!  end do
!  write(6, "('hprod_mklzmat. oay (I)')")
!  do ifun = 1, nfun
!     do jfun = 1, nfun
!        write(6, "(f20.10)", advance = 'no') aimag(oaz(jfun, ifun, 2))
!     end do
!     write(6, *)
!  end do
!  write(6, "('hprod_mklzmat. oaz (R)')")
!  do ifun = 1, nfun
!     do jfun = 1, nfun
!        write(6, "(f20.10)", advance = 'no') dble(oaz(jfun, ifun, 3))
!     end do
!     write(6, *)
!  end do
!  write(6, "('hprod_mklzmat. oaz (I)')")
!  do ifun = 1, nfun
!     do jfun = 1, nfun
!        write(6, "(f20.10)", advance = 'no') aimag(oaz(jfun, ifun, 3))
!     end do
!     write(6, *)
!  end do
  !DEBUG

  deallocate(oazp)

end subroutine hprod_mklzmat
!///////////////////////////////////////////////////////////////////////
