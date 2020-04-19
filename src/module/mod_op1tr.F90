!################################################################################
module mod_op1tr

  use, intrinsic :: iso_c_binding

  implicit none

  logical(c_bool), save :: do_op1tr = .false.

  complex(c_double_complex), save, allocatable :: op1tr_orb0(:,:,:)
  complex(c_double_complex), save, allocatable :: op1tr_orb(:,:,:)
  complex(c_double_complex), save, allocatable :: op1tr_dorb(:,:,:)
  complex(c_double_complex), save, allocatable :: op1tr_vorb(:,:,:)
  complex(c_double_complex), save, allocatable :: op1tr_aorb(:,:,:)
  complex(c_double_complex), save, allocatable :: op1tr_orbg(:,:,:)
  complex(c_double_complex), save, allocatable :: op1tr_dorbg(:,:,:)
  complex(c_double_complex), save, allocatable :: op1tr_vorbg(:,:,:)
  complex(c_double_complex), save, allocatable :: op1tr_aorbg(:,:,:)
  complex(c_double_complex), save, allocatable :: op1tr_umat(:,:)
  complex(c_double_complex), save, allocatable :: op1tr_dmat(:,:)
  complex(c_double_complex), save, allocatable :: op1tr_vmat(:,:)
  complex(c_double_complex), save, allocatable :: op1tr_amat(:,:)
  complex(c_double_complex), save, allocatable :: op1tr_dmat0(:,:)
  complex(c_double_complex), save, allocatable :: op1tr_vmat0(:,:)
  complex(c_double_complex), save, allocatable :: op1tr_amat0(:,:)
  complex(c_double_complex), save, allocatable :: den1_tr0(:,:)
  complex(c_double_complex), save, allocatable :: den1_tr1(:,:)
  complex(c_double_complex), save, allocatable :: den1_tr2(:,:)

  integer(c_int), save :: op1tr_nrad
  integer(c_int), pointer :: op1tr_nfun
  complex(c_double_complex), pointer :: op1tr_dP(:,:)
  complex(c_double_complex), pointer :: op1tr_vP(:,:)
  complex(c_double_complex), pointer :: op1tr_aP(:,:)
  complex(c_double_complex), pointer :: op1tr_dQ(:)
  complex(c_double_complex), pointer :: op1tr_vQ(:)
  complex(c_double_complex), pointer :: op1tr_aQ(:)

end module mod_op1tr
!################################################################################
