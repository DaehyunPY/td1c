!######################################################################
subroutine hprod_invden(den1, rden, rrden)

  use, intrinsic :: iso_c_binding
  use mod_control, only : throcc1, throcc2
  use mod_const, only : czero, runit, ctwo
  use mod_ormas, only : ncore, nact, nelact

  implicit none
  complex(c_double_complex), intent(in) :: den1(1:nact, 1:nact)
  complex(c_double_complex), intent(out) :: rden(1:nact, 1:nact)
  complex(c_double_complex), intent(out) :: rrden(1:nact, 1:nact)

  integer(c_int) :: iact, nact2
  complex(c_double_complex), allocatable :: tmp(:,:)

  if (nact == 0) return

  nact2 = nact * nact

!DEBUG  write(6, "('hprod_invden (0).')")
!DEBUG  do iact = 1, nact
!DEBUG     do jact = 1, nact
!DEBUG        write(6, "(f20.10)", advance = 'no') dble(den1(jact, iact))
!DEBUG     end do
!DEBUG     write(6, *)
!DEBUG  end do

  call hprod_invden_inv(nact, throcc1, den1, rden)

  if (ncore > 0) then
     if (nact * 2 == nelact(3)) then
        ! trivial active
        rrden(1:nact, 1:nact) = czero
     else
        ! non-trivial active
        allocate(tmp(nact, nact))
        call zcopy(nact2, czero, 0, tmp, 1)
        call zaxpy(nact2, -runit, den1, 1, tmp, 1)
        do iact = 1, nact
           tmp(iact, iact) = tmp(iact, iact) + ctwo
        end do
        call hprod_invden_inv(nact, throcc2, tmp, rrden)
        deallocate(tmp)
     end if
  end if

!debug
!  write(6, "('hprod_invden: rden')")
!!  write(6, "('hprod_invden: rrden')")
!  do iact = 1, nact
!     do jact = 1, nact
!        write(6, "(f20.10)", advance = 'no') abs(rden(jact, iact))
!!        write(6, "(f20.10)", advance = 'no') abs(rrden(jact, iact))
!     end do
!     write(6, *)
!  end do
!debug

end subroutine hprod_invden
!################################################################################
subroutine hprod_invden_inv(n, thresh, a, inva)

  use, intrinsic :: iso_c_binding
  use mod_control, only : reg_type
  use mod_const, only : czero, runit, zero, one

  implicit none
  integer(c_int), intent(in) :: n
  real(c_double), intent(in) :: thresh
  complex(c_double_complex), intent(in) :: a(1:n, 1:n)
  complex(c_double_complex), intent(out) :: inva(1:n, 1:n)

  integer(c_int) :: nsys
  integer(c_int) :: ifun, jfun, kfun
  complex(c_double_complex) :: diag, invd
  complex(c_double_complex), allocatable :: uvec(:,:)
  complex(c_double_complex), allocatable :: work(:,:)

  nsys = n
  allocate(uvec(1:n, 1:n))
  allocate(work(1:n, 1:n))

  call zcopy(nsys*nsys, a, 1, work, 1)
  call zcopy(nsys*nsys, czero, 0, inva, 1)
  call zcopy(nsys*nsys, czero, 0, uvec, 1)
  call lapack_zheev_dsc(nsys, work, uvec)
  
  do kfun = 1, n
     diag = work(kfun, kfun)
     if (thresh > zero) then
        invd = diag / (diag * diag + thresh)
     else
        invd = one / diag
     end if
     do jfun = 1, n
        do ifun = 1, n
           inva(ifun, jfun) = inva(ifun, jfun) + uvec(ifun, kfun) * invd &
                                       & * conjg(uvec(jfun, kfun))
        end do
     end do
  end do

  deallocate(work)
  deallocate(uvec)

end subroutine hprod_invden_inv
!######################################################################
