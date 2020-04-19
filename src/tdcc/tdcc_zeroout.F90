!######################################################################
subroutine tdcc_zeroout()

  use, intrinsic :: iso_c_binding
  use mod_cc, only : cc_rank,t0out,g0out,t1out,g1out,t2out,g2out,t3out,g3out
  implicit none

  t0out = 0d0
  if (cc_rank >= 1) t1out = 0d0
  if (cc_rank >= 2) t2out = 0d0
  if (cc_rank >= 3) t3out = 0d0

  g0out = 0d0
  if (cc_rank >= 1) g1out = 0d0
  if (cc_rank >= 2) g2out = 0d0
  if (cc_rank >= 3) g3out = 0d0

end subroutine tdcc_zeroout
!######################################################################
