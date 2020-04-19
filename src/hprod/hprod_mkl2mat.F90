!///////////////////////////////////////////////////////////////////////
subroutine hprod_mkl2mat(orb, oa2, oap, oam)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nfun
  use mod_sph, only : lmax1
  use mod_bas, only : mval
  use mod_rad, only : nrad, wrad
  use mod_const, only : one, half, czero, iunit

  implicit none
  complex(c_double_complex), intent(in) :: orb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(out) :: oa2(1:nfun, 1:nfun) ! l^2
  complex(c_double_complex), intent(out) :: oap(1:nfun, 1:nfun) ! l_+
  complex(c_double_complex), intent(out) :: oam(1:nfun, 1:nfun) ! l_-

  real(c_double) :: pcsi, pcsj, pcsji
  complex(c_double_complex) :: tmp1
  integer(c_long) :: ifun, jfun, irad, l, m, mi, mj, nproc, iproc, llr, ulr
  integer(c_long), external :: util_omp_iproc, util_omp_nproc
  complex(c_double_complex), allocatable :: oa2p(:,:,:) ! l^2
  complex(c_double_complex), allocatable :: oapp(:,:,:) ! l_+
  complex(c_double_complex), allocatable :: oamp(:,:,:) ! l_-

  nproc = util_omp_nproc()
  allocate(oa2p(1:nfun, 1:nfun, 0:(nproc-1)))
  allocate(oapp(1:nfun, 1:nfun, 0:(nproc-1)))
  allocate(oamp(1:nfun, 1:nfun, 0:(nproc-1)))
  oa2p(1:nfun, 1:nfun, 0:(nproc-1)) = czero
  oapp(1:nfun, 1:nfun, 0:(nproc-1)) = czero
  oamp(1:nfun, 1:nfun, 0:(nproc-1)) = czero

  !$omp parallel default(shared) private(iproc, llr, ulr, mi, mj, m, pcsi, pcsj, pcsji, tmp1)
  !###########################
  iproc = util_omp_iproc()
  call util_omp_disp(1, nrad - 1, llr, ulr)

  oa2p(1:nfun, 1:nfun, iproc) = czero
  oapp(1:nfun, 1:nfun, iproc) = czero
  oamp(1:nfun, 1:nfun, iproc) = czero

  do ifun = 1, nfun
     mi = mval(ifun)
     pcsi = (-one) ** max(0, mi)

     do jfun = 1, nfun
        mj = mval(jfun)
        pcsj = (-one) ** max(0, mj)
        pcsji = pcsj * pcsi
!       pcsji = one

        ! <lz>
        if (mj == mi) then
           m = mi
           do l = abs(m), lmax1
              tmp1 = czero
              do irad = llr, ulr
                 tmp1 = tmp1 + conjg(orb(irad, l, jfun)) * orb(irad, l, ifun)
              end do
              oa2p(jfun, ifun, iproc) = oa2p(jfun, ifun, iproc) + dble(l * (l + 1)) * tmp1
           end do
        end if

        ! <l+>
        if (mj == mi + 1) then
           m = mi
!          do l = abs(m), lmax1
           do l = abs(mj), lmax1
              tmp1 = czero
              do irad = llr, ulr
                 tmp1 = tmp1 + conjg(orb(irad, l, jfun)) * orb(irad, l, ifun)
              end do
              oapp(jfun, ifun, iproc) = oapp(jfun, ifun, iproc) + sqrt(dble((l - m)*(l + m + 1))) * tmp1
           end do
           oapp(jfun, ifun, iproc) = oapp(jfun, ifun, iproc) * pcsji
        end if

        ! <l->
        if (mj == mi - 1) then
           m = mi
!          do l = abs(m), lmax1
           do l = abs(mi), lmax1
              tmp1 = czero
              do irad = llr, ulr
                 tmp1 = tmp1 + conjg(orb(irad, l, jfun)) * orb(irad, l, ifun)
              end do
              oamp(jfun, ifun, iproc) = oamp(jfun, ifun, iproc) + sqrt(dble((l + m)*(l - m + 1))) * tmp1
           end do
           oamp(jfun, ifun, iproc) = oamp(jfun, ifun, iproc) * pcsji
        end if
     end do
  end do
  !###########################
  !$omp end parallel

  oa2(1:nfun, 1:nfun) = czero
  oap(1:nfun, 1:nfun) = czero
  oam(1:nfun, 1:nfun) = czero
  do iproc = 0, nproc - 1
     do ifun = 1, nfun
        do jfun = 1, nfun
           oa2(jfun, ifun) = oa2(jfun, ifun) + oa2p(jfun, ifun, iproc)
           oap(jfun, ifun) = oap(jfun, ifun) + oapp(jfun, ifun, iproc)
           oam(jfun, ifun) = oam(jfun, ifun) + oamp(jfun, ifun, iproc)
        end do
     end do
  end do

  !DEBUG
!  write(6, "('hprod_mkl2mat. oa2 (R)')")
!  do ifun = 1, nfun
!     do jfun = 1, nfun
!        write(6, "(f20.10)", advance = 'no') dble(oa2(jfun, ifun))
!     end do
!     write(6, *)
!  end do
!  write(6, "('hprod_mkl2mat. oa2 (I)')")
!  do ifun = 1, nfun
!     do jfun = 1, nfun
!        write(6, "(f20.10)", advance = 'no') aimag(oa2(jfun, ifun))
!     end do
!     write(6, *)
!  end do
!  write(6, "('hprod_mkl2mat. oa+ (R)')")
!  do ifun = 1, nfun
!     do jfun = 1, nfun
!        write(6, "(f20.10)", advance = 'no') dble(oap(jfun, ifun))
!     end do
!     write(6, *)
!  end do
!  write(6, "('hprod_mkl2mat. oa+ (I)')")
!  do ifun = 1, nfun
!     do jfun = 1, nfun
!        write(6, "(f20.10)", advance = 'no') aimag(oap(jfun, ifun))
!     end do
!     write(6, *)
!  end do
!  write(6, "('hprod_mkl2mat. oa- (R)')")
!  do ifun = 1, nfun
!     do jfun = 1, nfun
!        write(6, "(f20.10)", advance = 'no') dble(oam(jfun, ifun))
!     end do
!     write(6, *)
!  end do
!  write(6, "('hprod_mkl2mat. oa- (I)')")
!  do ifun = 1, nfun
!     do jfun = 1, nfun
!        write(6, "(f20.10)", advance = 'no') aimag(oam(jfun, ifun))
!     end do
!     write(6, *)
!  end do
  !DEBUG

  deallocate(oamp)
  deallocate(oapp)
  deallocate(oa2p)

end subroutine hprod_mkl2mat
!///////////////////////////////////////////////////////////////////////
