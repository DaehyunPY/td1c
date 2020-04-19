!################################################################################
subroutine ormas_arcwgt()

  use, intrinsic :: iso_c_binding

  use mod_ormas, only : iprint
  use mod_ormas, only : nelact, nact, arc_alph, arc_beta

  implicit none
  integer(c_long) :: m, k

  allocate(arc_alph(1:nact, 1:nelact(1)))
  allocate(arc_beta(1:nact, 1:nelact(2)))

  call ormas_arcwgt_spin(nelact(1), arc_alph)
  call ormas_arcwgt_spin(nelact(2), arc_beta)

  if (iprint > 2) then
     write(6,"('alpha arc weight:')")
     write(6, "(5x)", advance = 'no')
     do k = 1, nelact(1)
        write(6, "(i5)", advance = 'no') k
     end do
     write(6,*)
     do m = 1, nact
        write(6, "(i5)", advance = 'no') m
        do k = 1, nelact(1)
           write(6,"(i5)",advance='no') arc_alph(m, k)
        end do
        write(6,*)
     end do

     write(6,"('beta arc weight:')")
     write(6, "(5x)", advance = 'no')
     do k = 1, nelact(2)
        write(6, "(i5)", advance = 'no') k
     end do
     write(6,*)
     do m = 1, nact
        write(6, "(i5)", advance = 'no') m
        do k = 1, nelact(2)
           write(6,"(i5)",advance='no') arc_beta(m, k)
        end do
        write(6,*)
     end do
  end if

end subroutine ormas_arcwgt
!################################################################################
!################################################################################
subroutine ormas_arcwgt_spin(nel, arc)

  use, intrinsic :: iso_c_binding

  use mod_ormas, only : iprint
  use mod_ormas, only : nact

  implicit none
  integer(c_long), intent(in) :: nel
  integer(c_long), intent(inout) :: arc(1:nact, 1:*)
  integer(c_long), allocatable :: vtxwgt(:,:)
  integer(c_long) :: m, k

  ! initialization
  arc(1:nact, 1:nel) = 0

  allocate(vtxwgt(0:nact, 0:nel))
  vtxwgt(0:nact, 0:nel) = 0

  vtxwgt(0, 0) = 1
  do k = 0, nel
     vtxwgt(k, k) = 1
  end do     

  do m = 0, nact - nel
     vtxwgt(m, 0) = 1
  end do

  do k = 1, nel
     do m = k + 1, k + nact - nel
        vtxwgt(m, k) = vtxwgt(m-1, k) + vtxwgt(m-1, k-1)
     end do
  end do

  if (iprint > 2) then
     write(6,"('vertex weight:')")
     write(6, "(5x)", advance = 'no')
     do k = 0, nel
        write(6, "(i5)", advance = 'no') k
     end do
     write(6,*)
     do m = 0, nact
        write(6, "(i5)", advance = 'no') m
        do k = 0, nel
           write(6, "(i5)", advance = 'no') vtxwgt(m, k)
        end do
        write(6,*)
     end do
  end if

  do k = 1, nel
     do m = k, k + nact - nel
        arc(m, k) = vtxwgt(m, k) - vtxwgt(m-1, k-1)
     end do
  end do

  deallocate(vtxwgt)

end subroutine ormas_arcwgt_spin
!################################################################################
