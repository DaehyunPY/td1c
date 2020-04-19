!######################################################################
subroutine ormas_denipd(max_ipd,ovlp,cic,ipd,denipd)

  use,intrinsic :: iso_c_binding
  use mod_ormas,only : nact,nfun,ndetx,cic_old,tdcc

  implicit none
  !--------------------------------------------------------------------
  integer(c_int),intent(in) :: max_ipd
  complex(c_double_complex),intent(in) :: ovlp(1:nfun,1:nfun)
  complex(c_double_complex),intent(in) :: cic(1:ndetx)
  real(c_double),intent(out) :: ipd(0:max_ipd)
  complex(c_double),intent(out) :: denipd(1:nfun,1:nfun,0:max_ipd)
  !--------------------------------------------------------------------
  complex(c_double_complex), allocatable :: sout(:,:),D1(:,:),D2(:,:,:,:)
  integer(c_int) :: ifun

  allocate(sout(1:nfun,1:nfun))
  allocate(D1(1:nfun,1:nfun))
  allocate(D2(1:nfun,1:nfun,1:nfun,1:nfun))

  sout = -ovlp
  do ifun = 1, nfun
     sout(ifun, ifun) = sout(ifun, ifun) + 1d0
  end do
  call ormas_mkden1_full(cic,D1)
  call ormas_mkden2_full(cic,D2)

  denipd = 0d0
  call ormas_denipd_1(sout,D2,denipd(1,1,1))
  denipd(:,:,0) = D1(:,:) - denipd(:,:,1)

  deallocate(D2)
  deallocate(D1)
  deallocate(sout)

end subroutine ormas_denipd
!######################################################################
subroutine ormas_denipd_1(sout, D2, den1)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nfun

  implicit none
  complex(c_double_complex), intent(inout) :: sout(1:nfun, 1:nfun)
  complex(c_double_complex), intent(in) :: D2(1:nfun, 1:nfun, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(out) :: den1(1:nfun, 1:nfun)
  integer(c_int) :: ifun, jfun, kfun, lfun

  den1 = 0d0
  do ifun = 1, nfun
     do jfun = 1, nfun
        do kfun = 1, nfun
           do lfun = 1, nfun
              den1(ifun, jfun) = den1(ifun, jfun) + D2(ifun, jfun, kfun, lfun) * sout(lfun, kfun)
           end do
        end do
     end do
  end do

end subroutine ormas_denipd_1
!######################################################################
