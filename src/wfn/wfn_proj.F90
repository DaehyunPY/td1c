!///////////////////////////////////////////////////////////////////////
subroutine wfn_proj(porb, orb)

  use, intrinsic :: iso_c_binding

  use mod_bas, only : nbas
  use mod_ormas, only : nocc, nvir, nfun

  implicit none
  complex(c_double_complex), intent(in) :: porb(1:nbas, 1:nfun)
  complex(c_double_complex), intent(inout) :: orb(1:nbas, 1:nfun)

  complex(c_double_complex) :: fac, fac1, fac2
  complex(c_double_complex), allocatable :: s2p(:,:)
  complex(c_double_complex), allocatable :: s2a(:,:)
  integer(c_int) :: ifun, jfun, ibas

  write(6, "('wfn_proj may include bug!')")
  stop

  allocate(s2p(1:nfun, 1:nfun))
  if (nvir > 0) allocate(s2a(1:nfun, 1:nfun))

  ! <p|F_i|i> {p} \in {i,a}
  do ifun = 1, nfun
     do jfun = 1, nfun
        fac = (0.d+0, 0.d+0)
        do ibas = 1, nbas
           fac = fac + conjg(porb(ibas, jfun)) * orb(ibas, ifun)
        end do
        s2p(jfun, ifun) = fac
     end do
  end do

  ! <a|F_i|i>
  if (nvir > 0) then
     do ifun = 1, nocc
        do jfun = nocc + 1, nfun
           fac = (0.d+0, 0.d+0)
           do ibas = 1, nbas
              fac = fac + conjg(porb(ibas, jfun)) * orb(ibas, ifun)
           end do
           s2a(jfun, ifun) = fac
        end do
     end do
  end if

!debug  do ifun = 1, nfun
!debug     do jfun = 1, nfun
!debug        write(6, "('wfn_proj: s2p', 2i5, 2f20.10)") jfun, ifun, s2p(jfun, ifun)
!debug     end do
!debug  end do

  ! 1 - \sum_p |p><p|F_i|i>, {p} \in {i,a}
  do ifun = 1, nfun
     do jfun = 1, nfun
        fac = s2p(jfun, ifun)
        do ibas = 1, nbas
           orb(ibas, ifun) = orb(ibas, ifun) - porb(ibas, jfun) * fac
        end do
     end do
  end do

  ! + \sum_a |a><a|F_i|i>
  if (nvir > 0) then
     do ifun = 1, nocc
        do jfun = nocc + 1, nfun
           fac1 = s2a(jfun, ifun)
           fac2 = conjg(fac1)
           do ibas = 1, nbas
              orb(ibas, ifun) = orb(ibas, ifun) + porb(ibas, jfun) * fac1
              orb(ibas, jfun) = orb(ibas, jfun) - porb(ibas, ifun) * fac2
           end do
        end do
     end do
  end if

  if (nvir > 0) deallocate(s2a)
  deallocate(s2p)

end subroutine wfn_proj
!///////////////////////////////////////////////////////////////////////
