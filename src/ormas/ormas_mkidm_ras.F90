!################################################################################
subroutine ormas_mkidm_ras(ovlp, cic, dens)

  use, intrinsic :: iso_c_binding
  use mod_const, only : zero, one, czero, ctwo
  use mod_ormas, only : thrdet, nelact, neltot, ncore, nact, nfun, nstr_alph, nstr_beta, &
       & orb_alph, orb_beta, onv_alph, onv_beta

  implicit none
  !--------------------------------------------------------------------
  complex(c_double_complex), intent(in) :: ovlp(1:nfun, 1:nfun)
  complex(c_double_complex), intent(in) :: cic(1:nstr_alph, 1:nstr_beta)
  complex(c_double_complex), intent(out) :: dens(1:nfun, 1:nfun)
  !--------------------------------------------------------------------
  complex(c_double_complex) :: dets
  complex(c_double_complex), external :: ormas_dets0, ormas_dets_single
  complex(c_double_complex), allocatable :: work(:,:)
  complex(c_double_complex), allocatable :: work1(:,:)
  complex(c_double_complex), allocatable :: cdeta(:,:)
  complex(c_double_complex), allocatable :: detbc(:,:)
  complex(c_double_complex), allocatable :: cdeta1(:,:)
  complex(c_double_complex), allocatable :: detbc1(:,:)
  complex(c_double_complex), allocatable :: denp(:,:,:)
  integer(c_long), external :: util_omp_nproc
  integer(c_long), external :: util_omp_iproc
  integer(c_long) :: ndet, istr, jstr, ifun, jfun, ii, iproc, nproc, maxne2

  nproc = util_omp_nproc()
  ndet = nstr_alph * nstr_beta
  maxne2 = max(neltot(1), neltot(2)) ** 2
  allocate(work(1:maxne2, 0:nproc-1))
  allocate(work1(1:maxne2, 0:nproc-1))
  allocate(cdeta(1:nstr_beta, 1:nstr_alph))
  allocate(detbc(1:nstr_beta, 1:nstr_alph))
  allocate(cdeta1(1:nstr_beta, 0:nproc-1))
  allocate(detbc1(1:nstr_alph, 0:nproc-1))
  allocate(denp(1:nfun, 1:nfun, 0:nproc-1))

  cdeta(1:nstr_beta, 1:nstr_alph) = czero
  detbc(1:nstr_beta, 1:nstr_alph) = czero

  !$omp parallel default(shared) private(istr,jstr,dets,ii,iproc)
  iproc = util_omp_iproc()
  !$omp do
  do jstr = 1, nstr_alph
     do istr = 1, nstr_alph
        dets = ormas_dets0(istr, jstr, ncore, neltot(1), nelact(1), nfun, orb_alph, ovlp, work(1,iproc))
        if (abs(dets) > thrdet) then
           do ii = 1, nstr_beta
              cdeta(ii, jstr) = cdeta(ii, jstr) + conjg(cic(istr, ii)) * dets
           end do
        end if
     end do
  end do
  !$omp end do
  !$omp do
  do istr = 1, nstr_beta
     do jstr = 1, nstr_beta
        dets = ormas_dets0(istr, jstr, ncore, neltot(2), nelact(2), nfun, orb_beta, ovlp, work(1,iproc))
        if (abs(dets) > thrdet) then
           do ii = 1, nstr_alph
              detbc(istr, ii) = detbc(istr, ii) + dets * cic(ii, jstr)
           end do
        end if
     end do
  end do
  !$omp end do
  !$omp end parallel

  !$omp parallel default(shared) private(istr,jstr,dets,ii,iproc)
  iproc = util_omp_iproc()
  denp(1:nfun, 1:nfun, iproc) = czero
  !$omp do
  do jstr = 1, nstr_alph
     do ifun = 1, nfun
        do jfun = 1, nfun
           cdeta1(1:nstr_beta, iproc) = czero
           do istr = 1, nstr_alph
              dets = ormas_dets0(istr, jstr, ncore, neltot(1), nelact(1), nfun, orb_alph, ovlp, work(1,iproc))
              dets = ormas_dets_single(ifun, jfun, neltot(1), onv_alph(1,istr), onv_alph(1,jstr), work(1,iproc), work1(1,iproc))
              if (abs(dets) > thrdet) then
                 do ii = 1, nstr_beta
                    cdeta1(ii, iproc) = cdeta1(ii, iproc) + conjg(cic(istr, ii)) * dets
                 end do
              end if
           end do
           do ii = 1, nstr_beta
              denp(ifun, jfun, iproc) = denp(ifun, jfun, iproc) + cdeta1(ii, iproc) * detbc(ii, jstr)
           end do
        end do
     end do
  end do
  !$omp end do

  !$omp do
  do istr = 1, nstr_beta
     do ifun = 1, nfun
        do jfun = 1, nfun
           detbc1(1:nstr_alph, iproc) = czero
           do jstr = 1, nstr_beta
              dets = ormas_dets0(istr, jstr, ncore, neltot(2), nelact(2), nfun, orb_beta, ovlp, work(1,iproc))
              dets = ormas_dets_single(ifun, jfun, neltot(2), onv_beta(1,istr), onv_beta(1,jstr), ovlp, work(1,iproc))
              if (abs(dets) > thrdet) then
                 do ii = 1, nstr_alph
                    detbc1(ii, iproc) = detbc1(ii, iproc) + dets * cic(ii, jstr)
                 end do
              end if
           end do
           do ii = 1, nstr_alph
              denp(ifun, jfun, iproc) = denp(ifun, jfun, iproc) + cdeta(jstr, ii) * detbc1(ii, iproc)
           end do
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel

  dens(1:nfun, 1:nfun) = czero
  do iproc = 0, nproc - 1
     dens(1:nfun, 1:nfun) = dens(1:nfun, 1:nfun) + denp(1:nfun, 1:nfun, iproc)
  end do

  !debug
!  write(6, "('ormas_mkidm_ras: Ovlp-R')")
!  do ifun = 1, nfun
!     do jfun = 1, nfun
!        write(6, "(f15.8)", advance = 'no') dble(ovlp(jfun, ifun))
!     end do
!     write(6, *)
!  end do
!  write(6, "('ormas_mkidm_ras: Ovlp-I')")
!  do ifun = 1, nfun
!     do jfun = 1, nfun
!        write(6, "(f15.8)", advance = 'no') aimag(ovlp(jfun, ifun))
!     end do
!     write(6, *)
!  end do
!  write(6, "('ormas_mkidm_ras: IDM-R')")
!  do ifun = 1, nfun
!     do jfun = 1, nfun
!        write(6, "(f15.8)", advance = 'no') dble(dens(jfun, ifun))
!     end do
!     write(6, *)
!  end do
!  write(6, "('ormas_mkidm_ras: IDM-I')")
!  do ifun = 1, nfun
!     do jfun = 1, nfun
!        write(6, "(f15.8)", advance = 'no') aimag(dens(jfun, ifun))
!     end do
!     write(6, *)
!  end do
!  stop
  !debug

  deallocate(denp)
  deallocate(detbc1)
  deallocate(cdeta1)
  deallocate(detbc)
  deallocate(cdeta)
  deallocate(work1)
  deallocate(work)

end subroutine ormas_mkidm_ras
!################################################################################
