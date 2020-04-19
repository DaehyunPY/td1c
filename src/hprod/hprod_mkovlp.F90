!///////////////////////////////////////////////////////////////////////
subroutine hprod_mkovlp(rmax, orb1, orb2, ovlp)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad, xrad
  use mod_sph, only : lmax1
  use mod_ormas, only : nfun
  use mod_bas, only : mval
  use mod_const, only : czero

  implicit none
  real(c_double), intent(in) :: rmax
  complex(c_double_complex), intent(in) :: orb1(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: orb2(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(out) :: ovlp(1:nfun, 1:nfun)

  integer(c_long) :: ifun, jfun, max_irad, irad, l, mi, mj, m, nproc, iproc, llr, ulr
  integer(c_long), external :: util_omp_nproc, util_omp_iproc
  complex(c_double_complex) :: tmp
  complex(c_double_complex), allocatable :: ovlpp(:,:,:)

  nproc = util_omp_nproc()
  allocate(ovlpp(1:nfun, 1:nfun, 0:(nproc-1)))
  ovlp(1:nfun, 1:nfun) = czero
  ovlpp(1:nfun, 1:nfun, 0:(nproc-1)) = czero

  ! determine max_irad
  max_irad = nrad - 1
  do irad = 1, nrad - 1
     if (xrad(irad) > rmax) then
        max_irad = irad
        exit
     end if
  end do

  !$omp parallel default(shared) private(iproc, mi, mj, m, tmp, llr, ulr)
  !###########################
  iproc = util_omp_iproc()
  ovlpp(1:nfun, 1:nfun, iproc) = czero

  call util_omp_disp(1, max_irad, llr, ulr)
  do ifun = 1, nfun
     mi = mval(ifun)
     do jfun = 1, nfun
        mj = mval(jfun)
        if (mi == mj) then
           m = mi
           tmp = czero
           do l = abs(m), lmax1
              do irad = llr, ulr
                 tmp = tmp + conjg(orb1(irad, l, jfun)) &
                                 * orb2(irad, l, ifun)
              end do
           end do
           ovlpp(jfun, ifun, iproc) = tmp
        end if
     end do
  end do
  !###########################
  !$omp end parallel

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

  !DEBUG
  !write(6, "('hprod_mkovlp. (R)')")
  !do ifun = 1, nfun
  !   do jfun = 1, nfun
  !      write(6, "(f20.10)") dble(ovlp(jfun, ifun))
  !   end do
  !   write(6, *)
  !end do
  !write(6, "('hprod_mkovlp. (I)')")
  !do ifun = 1, nfun
  !   do jfun = 1, nfun
  !      write(6, "(f20.10)") aimag(ovlp(jfun, ifun))
  !   end do
  !   write(6, *)
  !end do
  !DEBUG

end subroutine hprod_mkovlp
!///////////////////////////////////////////////////////////////////////
subroutine hprod_mkovlp_shell(rmin, rmax, orb1, orb2, ovlp)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad, xrad
  use mod_sph, only : lmax1
  use mod_ormas, only : nfun
  use mod_bas, only : mval
  use mod_const, only : czero

  implicit none
  real(c_double), intent(in) :: rmin, rmax
  complex(c_double_complex), intent(in) :: orb1(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: orb2(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(out) :: ovlp(1:nfun, 1:nfun)

  integer(c_long) :: ifun, jfun, min_irad, max_irad, irad, l, mi, mj, m, nproc, iproc, llr, ulr
  integer(c_long), external :: util_omp_nproc, util_omp_iproc
  complex(c_double_complex) :: tmp
  complex(c_double_complex), allocatable :: ovlpp(:,:,:)

  nproc = util_omp_nproc()
  allocate(ovlpp(1:nfun, 1:nfun, 0:(nproc-1)))
  ovlp(1:nfun, 1:nfun) = czero
  ovlpp(1:nfun, 1:nfun, 0:(nproc-1)) = czero

  ! determine min_irad
  min_irad = 1
  do irad = nrad - 1, 1, -1
     if (xrad(irad) < rmin) then
        min_irad = irad
        exit
     end if
  end do

  ! determine max_irad
  max_irad = nrad - 1
  do irad = 1, nrad - 1
     if (xrad(irad) > rmax) then
        max_irad = irad
        exit
     end if
  end do

  if (min_irad > max_irad) stop 'hprod_mkovlp_shell: min_irad > max_irad.'

  !$omp parallel default(shared) private(iproc, mi, mj, m, tmp, llr, ulr)
  !###########################
  iproc = util_omp_iproc()
  ovlpp(1:nfun, 1:nfun, iproc) = czero

  call util_omp_disp(min_irad, max_irad, llr, ulr)
  do ifun = 1, nfun
     mi = mval(ifun)
     do jfun = 1, nfun
        mj = mval(jfun)
        if (mi == mj) then
           m = mi
           tmp = czero
           do l = abs(m), lmax1
              do irad = llr, ulr
                 tmp = tmp + conjg(orb1(irad, l, jfun)) &
                                 * orb2(irad, l, ifun)
              end do
           end do
           ovlpp(jfun, ifun, iproc) = tmp
        end if
     end do
  end do
  !###########################
  !$omp end parallel

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

  !DEBUG
  ! write(6, "('hprod_mkovlp. (R)')")
  ! do ifun = 1, nfun
  !    do jfun = 1, nfun
  !       write(6, "(f20.10)") dble(ovlp(jfun, ifun))
  !    end do
  !    write(6, *)
  ! end do
  ! write(6, "('hprod_mkovlp. (I)')")
  ! do ifun = 1, nfun
  !    do jfun = 1, nfun
  !       write(6, "(f20.10)") aimag(ovlp(jfun, ifun))
  !    end do
  !    write(6, *)
  ! end do
  !DEBUG

end subroutine hprod_mkovlp_shell
!///////////////////////////////////////////////////////////////////////
subroutine hprod_mkovlp_kshell(kmin, kmax, orb1, orb2, ovlp)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad, xrad
  use mod_sph, only : lmax1
  use mod_ormas, only : nfun
  use mod_bas, only : mval
  use mod_const, only : czero
  use mod_pes, only : pes_numk, pes_k_min, pes_k_step

  implicit none
  real(c_double), intent(in) :: kmin, kmax
  complex(c_double_complex), intent(in) :: orb1(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: orb2(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(out) :: ovlp(1:nfun, 1:nfun)

  integer(c_long) :: ifun, jfun, min_ik, max_ik, ik, l, mi, mj, m, nproc, iproc, llk, ulk
  integer(c_long), external :: util_omp_nproc, util_omp_iproc
  real(c_double) :: kval
  complex(c_double_complex) :: tmp
  complex(c_double_complex), allocatable :: ovlpp(:,:,:)
  complex(c_double_complex), allocatable :: orbk1(:,:,:), orbk2(:,:,:)


  nproc = util_omp_nproc()
  allocate(orbk1(0:pes_numk, 0:lmax1, 1:nfun))
  allocate(orbk2(0:pes_numk, 0:lmax1, 1:nfun))
  allocate(ovlpp(1:nfun, 1:nfun, 0:(nproc-1)))
  ovlp(1:nfun, 1:nfun) = czero
  ovlpp(1:nfun, 1:nfun, 0:(nproc-1)) = czero

  call pes_spec1_k_psik(orb1, orbk1)
  call pes_spec1_k_psik(orb2, orbk2)


  ! determine min_ik
  min_ik = 0
  do ik = pes_numk, 0, -1
     kval = pes_k_min + ik * pes_k_step
     if (kval < kmin) then
        min_ik = ik
        exit
     end if
  end do

  ! determine max_ik
  max_ik = pes_numk
  do ik = 0, pes_numk
     kval = pes_k_min + ik * pes_k_step
     if (kval > kmax) then
        max_ik = ik
        exit
     end if
  end do

  !$omp parallel default(shared) private(iproc,mi,mj,m,tmp,llk,ulk)
  !###########################
  iproc = util_omp_iproc()
  call util_omp_disp(min_ik, max_ik, llk, ulk)
  do ifun = 1, nfun
     mi = mval(ifun)
     do jfun = 1, nfun
        mj = mval(jfun)
        if (mi == mj) then
           m = mi
           tmp = czero
           do l = abs(m), lmax1
              do ik = llk, ulk
                 tmp = tmp + conjg(orbk1(ik,l,ifun))*orbk2(ik,l,jfun)
              end do
           end do
           ovlpp(ifun,jfun,iproc) = ovlpp(ifun,jfun,iproc) + tmp
        end if
     end do
  end do
  !###########################
  !$omp end parallel

  do iproc = 0, nproc - 1
     do ifun = 1, nfun
        do jfun = 1, nfun
           ovlp (jfun, ifun) = &
         & ovlp (jfun, ifun) + &
         & ovlpp(jfun, ifun, iproc)
        end do
     end do
  end do

  deallocate(orbk1)
  deallocate(orbk2)
  deallocate(ovlpp)

  !DEBUG
  ! write(6, "('hprod_mkovlp. (R)')")
  ! do ifun = 1, nfun
  !    do jfun = 1, nfun
  !       write(6, "(f20.10)") dble(ovlp(jfun, ifun))
  !    end do
  !    write(6, *)
  ! end do
  ! write(6, "('hprod_mkovlp. (I)')")
  ! do ifun = 1, nfun
  !    do jfun = 1, nfun
  !       write(6, "(f20.10)") aimag(ovlp(jfun, ifun))
  !    end do
  !    write(6, *)
  ! end do
  !DEBUG

end subroutine hprod_mkovlp_kshell
!///////////////////////////////////////////////////////////////////////
