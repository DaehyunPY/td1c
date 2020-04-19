!################################################################################
module mod_cc

  use, intrinsic :: iso_c_binding
  implicit none

  public

  integer(c_long), save, pointer :: cc_type
  character(len=16) :: cc_code
  integer(c_long) :: cc_rank, norb1, norb2
  logical(c_bool) :: optcc, bcc, optbcc
  logical(c_bool) :: dot1

  logical(c_bool), parameter :: tonly = .false.
  integer(c_long), parameter :: cc_maxcyc = 40
  logical(c_bool), parameter :: cc_nonredundant = .true.

  logical(c_bool), parameter :: cc_solve = .true.
  logical(c_bool), parameter :: cc_solve_itr = .false.
  logical(c_bool), parameter :: cc_sreal = .true.
  logical(c_bool), parameter :: cc_read_ort = .false.
  integer(c_long), parameter :: cc_xij_xab = 0
  integer(c_long), parameter :: nspin = 1
  integer(c_long), parameter :: nbiort = 1
  real(c_double), parameter :: thrtamp = 1D-15
  real(c_double), parameter :: thrgamp = 1D-15

  integer(c_long) :: ncc0,ncc1a,ncc2aa,ncc2ab,ncc3aaa,ncc3aab
  integer(c_long),allocatable :: map_cc0(:,:)
  integer(c_long),allocatable :: map_cc1a(:,:)
  integer(c_long),allocatable :: map_cc2aa(:,:)
  integer(c_long),allocatable :: map_cc2ab(:,:)
  integer(c_long),allocatable :: map_cc3aaa(:,:)
  integer(c_long),allocatable :: map_cc3aab(:,:)
  integer(c_long),allocatable :: h1_cc1a(:),p1_cc1a(:)
  integer(c_long),allocatable :: h1_cc2aa(:),h2_cc2aa(:),p1_cc2aa(:),p2_cc2aa(:)
  integer(c_long),allocatable :: h1_cc2ab(:),h2_cc2ab(:),p1_cc2ab(:),p2_cc2ab(:)
  integer(c_long),allocatable :: h1_cc3aaa(:),h2_cc3aaa(:),h3_cc3aaa(:),p1_cc3aaa(:),p2_cc3aaa(:),p3_cc3aaa(:)
  integer(c_long),allocatable :: h1_cc3aab(:),h2_cc3aab(:),h3_cc3aab(:),p1_cc3aab(:),p2_cc3aab(:),p3_cc3aab(:)

  complex(c_double_complex), allocatable :: cc_i0(:)
  complex(c_double_complex), allocatable :: cc_i1(:)
  complex(c_double_complex), allocatable :: cc_i2(:)
  complex(c_double_complex), allocatable :: cc_i3(:)

  complex(c_double_complex), allocatable :: fock(:,:,:)
  complex(c_double_complex), allocatable :: int2x(:,:,:,:,:)

  complex(c_double_complex) :: t0inp, dt0inp, t0out, dt0out
  complex(c_double_complex) :: g0inp, dg0inp, g0out, dg0out
  complex(c_double_complex), allocatable :: t1inp(:,:,:), t1out(:,:,:)
  complex(c_double_complex), allocatable :: g1inp(:,:,:), g1out(:,:,:)
  complex(c_double_complex), allocatable :: t2inp(:,:,:,:,:), t2out(:,:,:,:,:)
  complex(c_double_complex), allocatable :: g2inp(:,:,:,:,:), g2out(:,:,:,:,:)
  complex(c_double_complex), allocatable :: t3inp(:,:,:,:,:,:,:), t3out(:,:,:,:,:,:,:)
  complex(c_double_complex), allocatable :: g3inp(:,:,:,:,:,:,:), g3out(:,:,:,:,:,:,:)
  complex(c_double_complex), allocatable :: dt1inp(:,:,:), dt1out(:,:,:)
  complex(c_double_complex), allocatable :: dg1inp(:,:,:), dg1out(:,:,:)
  complex(c_double_complex), allocatable :: dt2inp(:,:,:,:,:), dt2out(:,:,:,:,:)
  complex(c_double_complex), allocatable :: dg2inp(:,:,:,:,:), dg2out(:,:,:,:,:)
  complex(c_double_complex), allocatable :: dt3inp(:,:,:,:,:,:,:), dt3out(:,:,:,:,:,:,:)
  complex(c_double_complex), allocatable :: dg3inp(:,:,:,:,:,:,:), dg3out(:,:,:,:,:,:,:)

  complex(c_double_complex), allocatable :: den1s(:,:,:)
  complex(c_double_complex), allocatable :: den2s(:,:,:,:,:)
  complex(c_double_complex), allocatable :: den1_noref(:,:,:)
  complex(c_double_complex), allocatable :: den2_noref(:,:,:,:,:)

  complex(c_double_complex), allocatable :: itm_hh(:,:)
  complex(c_double_complex), allocatable :: itm_pp(:,:)
  complex(c_double_complex), allocatable :: itm_hp(:,:)
  complex(c_double_complex), allocatable :: itm_ph(:,:)
  complex(c_double_complex), allocatable :: itm_hhhh(:,:,:,:)
  complex(c_double_complex), allocatable :: itm_pppp(:,:,:,:)
  complex(c_double_complex), allocatable :: itm_pphp(:,:,:,:)
  complex(c_double_complex), allocatable :: itm_hphp(:,:,:,:)
  complex(c_double_complex), allocatable :: itm_hhhp(:,:,:,:)
  complex(c_double_complex), allocatable :: itm_hphh(:,:,:,:)
  complex(c_double_complex), allocatable :: itm_hpph(:,:,:,:)

  complex(c_double_complex), allocatable :: itm_hhh(:,:,:)
  complex(c_double_complex), allocatable :: itm_hph(:,:,:)
  complex(c_double_complex), allocatable :: itm_pph(:,:,:)
  complex(c_double_complex), allocatable :: itm_ppp(:,:,:)
  complex(c_double_complex), allocatable :: itm_hppp(:,:,:,:)
  complex(c_double_complex), allocatable :: itm_hhph(:,:,:,:)
  complex(c_double_complex), allocatable :: itm_phpp(:,:,:,:)
  complex(c_double_complex), allocatable :: itm_ppph(:,:,:,:)
  complex(c_double_complex), allocatable :: itm_hhphhh(:,:,:,:,:,:)

end module mod_cc
!################################################################################
