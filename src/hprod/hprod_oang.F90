!///////////////////////////////////////////////////////////////////////
subroutine hprod_oang(orb, cic, oangz, oang2)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : mval
  use mod_const, only : zero
  use mod_ormas, only : ncore, nact
  use mod_hprod, only : den1, den2

  implicit none
  complex(c_double_complex), intent(in) :: orb(1:*)
  complex(c_double_complex), intent(in) :: cic(1:*)
  real(c_double), intent(out) :: oangz, oang2

  call ormas_mkden1(cic, den1);
  call ormas_mkden2(cic, den1, den2);
  call hprod_oangz(den1, oangz);
  call hprod_oang2(orb, den1, den2, oang2);

end subroutine hprod_oang
!///////////////////////////////////////////////////////////////////////
subroutine hprod_oangz(den1, oangz)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : ncore, nact
  use mod_bas, only : mval
  use mod_const, only : zero

  implicit none
  !##### arg #####
  complex(c_double_complex), intent(in) :: den1(1:nact, 1:nact)
  real(c_double), intent(out) :: oangz
  !##### local #####
  integer(c_long) :: iact, ifun

  oangz = zero
  do ifun = 1, ncore
     oangz = oangz + dble(mval(ifun))
  end do
  oangz = oangz + oangz

  do iact = 1, nact
     ifun = ncore + iact
     oangz = oangz + dble(mval(ifun) * den1(iact, iact))
  end do

end subroutine hprod_oangz
!///////////////////////////////////////////////////////////////////////
subroutine hprod_oang2(orb, den1, den2, oang2)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : mval
  use mod_ormas, only : ncore, nact, nfun
  use mod_const, only : one, half, zero, two
  use mod_hprod, only : int1e, int2e

  implicit none
  !##### args #####
  complex(c_double_complex), intent(in) :: orb(1:*)
  complex(c_double_complex), intent(in) :: den1(1:nact, 1:nact)
  complex(c_double_complex), intent(in) :: den2(1:nact, 1:nact, 1:nact, 1:nact)
  real(c_double), intent(out) :: oang2
  !##### local #####
  integer(c_long) :: ifun, jfun, kfun, lfun
  integer(c_long) :: iact, jact, kact, lact
  real(c_double) :: m_i, m_k, oang2_core, oang2_act1, oang2_act2
  complex(c_double_complex), allocatable :: oa2(:,:) ! l^2
  complex(c_double_complex), allocatable :: oap(:,:) ! l_+
  complex(c_double_complex), allocatable :: oam(:,:) ! l_-

  allocate(oa2(1:nfun, 1:nfun))
  allocate(oap(1:nfun, 1:nfun))
  allocate(oam(1:nfun, 1:nfun))
  call hprod_mkl2mat(orb, oa2, oap, oam)

  !NYI
  oang2_core = zero
  !NYI

  do iact = 1, nact
     ifun = ncore + iact
     do jact = 1, nact
        jfun = ncore + jact
        int1e(jact, iact) = oa2(jfun, ifun)
     end do
  end do
  do iact = 1, nact
     ifun = ncore + iact
     do jact = 1, nact
        jfun = ncore + jact
        do kact = 1, nact
           kfun = ncore + kact
           do lact = 1, nact
              lfun = ncore + lact
              int2e(lact, kact, jact, iact) = oap(lfun, kfun) * oam(jfun, ifun) &
                                            + oam(lfun, kfun) * oap(jfun, ifun)
           end do
        end do
     end do
  end do
  do iact = 1, nact
     ifun = ncore + iact
     m_i = dble(mval(ifun))
     do kact = 1, nact
        kfun = ncore + kact
        m_k = dble(mval(kfun))
        int2e(kact, kact, iact, iact) = int2e(kact, kact, iact, iact) + m_i * m_k * two
     end do
  end do

  oang2_act1 = zero
  oang2_act2 = zero
  call hprod_op1e(one, int1e, den1, oang2_act1)
  call hprod_op2e(half, int2e, den2, oang2_act2)
  !DEBUG
  !DEBUG write(6, "('hprod_oang2: L^2_C    = ', f20.10)") oang2_core
  !DEBUG write(6, "('hprod_oang2: L^2_A(O) = ', f20.10)") oang2_act1
  !DEBUG write(6, "('hprod_oang2: L^2_A(T) = ', f20.10)") oang2_act2
  !DEBUG

  oang2 = oang2_core + oang2_act1 + oang2_act2

  deallocate(oam)
  deallocate(oap)
  deallocate(oa2)

end subroutine hprod_oang2
!///////////////////////////////////////////////////////////////////////
