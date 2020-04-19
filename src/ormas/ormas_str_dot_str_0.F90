!################################################################################
subroutine ormas_str_dot_str_0(iipx,istr,jstr,ncore,nel,nela,orb,s0,work,ovlp)

  use,intrinsic :: iso_c_binding
  use mod_ormas,only : thrdet,nfun
  use mod_const,only : czero,runit

  implicit none
  !--------------------------------------------------------------------
  integer(c_int),intent(in) :: iipx,istr,jstr,ncore,nel,nela,orb(0:nela,*)
  complex(c_double_complex),intent(in) :: s0(1:nel,1:nel)
  complex(c_double_complex),intent(out) :: work(1:(nel-iipx),1:(nel-iipx))
  complex(c_double_complex),intent(out) :: ovlp
  !--------------------------------------------------------------------

  if (iipx == 0) then
     call ormas_str_dot_str_0_0(istr,jstr,ncore,nel,nela,orb,s0,work,ovlp)
  else if (iipx == 1) then
     call ormas_str_dot_str_0_1(istr,jstr,ncore,nel,nela,orb,s0,work,ovlp)
  else if (iipx == 2) then
     call ormas_str_dot_str_0_2(istr,jstr,ncore,nel,nela,orb,s0,work,ovlp)
  else if (iipx == 3) then
     call ormas_str_dot_str_0_3(istr,jstr,ncore,nel,nela,orb,s0,work,ovlp)
  else if (iipx == 4) then
     call ormas_str_dot_str_0_4(istr,jstr,ncore,nel,nela,orb,s0,work,ovlp)
  else
     stop 'ormas_str_dot_str_0: bad iipx.'
  end if

  return

end subroutine ormas_str_dot_str_0
!################################################################################
subroutine ormas_str_dot_str_0_0(istr,jstr,ncore,nel,nela,orb,s0,work,ovlp)

  use,intrinsic :: iso_c_binding
  use mod_ormas,only : thrdet,nfun
  use mod_const,only : czero,runit

  implicit none
  !--------------------------------------------------------------------
  integer(c_int),intent(in) :: istr,jstr,ncore,nel,nela,orb(0:nela,*)
  complex(c_double_complex),intent(in) :: s0(1:nel,1:nel)
  complex(c_double_complex),intent(out) :: work(1:nel,1:nel)
  complex(c_double_complex),intent(out) :: ovlp
  !--------------------------------------------------------------------
  complex(c_double_complex),external :: util_det

  ovlp = czero
  work = s0
  ovlp = ovlp + util_det(nel,thrdet,work)

end subroutine ormas_str_dot_str_0_0
!################################################################################
subroutine ormas_str_dot_str_0_1(istr,jstr,ncore,nel,nela,orb,s0,work,ovlp)

  use,intrinsic :: iso_c_binding
  use mod_ormas,only : thrdet,nfun
  use mod_const,only : czero,runit

  implicit none
  !--------------------------------------------------------------------
  integer(c_int),intent(in) :: istr,jstr,ncore,nel,nela,orb(0:nela,*)
  complex(c_double_complex),intent(in) :: s0(1:nel,1:nel)
  complex(c_double_complex),intent(out) :: work(1:(nel-1),1:(nel-1))
  complex(c_double_complex),intent(out) :: ovlp
  !--------------------------------------------------------------------
  integer(c_int) :: i1,j1,sgn1
  complex(c_double_complex),external :: util_det
  logical(c_bool),external :: ormas_chkorth

  ovlp = czero
  do i1 = 1,nel
     do j1 = 1,nel
        if (.not.ormas_chkorth(i1,j1,istr,jstr,nela,orb)) cycle
        sgn1 = (-1)**(i1+j1)
        call util_window1(nel,i1,j1,s0,work)
        ovlp = ovlp + util_det(nel-1,thrdet,work)*sgn1
     end do
  end do

end subroutine ormas_str_dot_str_0_1
!################################################################################
subroutine ormas_str_dot_str_0_2(istr,jstr,ncore,nel,nela,orb,s0,work,ovlp)

  use,intrinsic :: iso_c_binding
  use mod_ormas,only : thrdet,nfun
  use mod_const,only : czero,runit

  implicit none
  !--------------------------------------------------------------------
  integer(c_int),intent(in) :: istr,jstr,ncore,nel,nela,orb(0:nela,*)
  complex(c_double_complex),intent(in) :: s0(1:nel,1:nel)
  complex(c_double_complex),intent(out) :: work(1:(nel-2),1:(nel-2))
  complex(c_double_complex),intent(out) :: ovlp
  !--------------------------------------------------------------------
  integer(c_int) :: i1,i2,j1,j2,sgn1,sgn12
  complex(c_double_complex),external :: util_det
  logical(c_bool),external :: ormas_chkorth

  ovlp = czero
  do i1 = 1,nel
     do j1 = 1,nel
        if (.not.ormas_chkorth(i1,j1,istr,jstr,nela,orb)) cycle
        sgn1 = (-1)**(i1+j1)
        do i2 = i1 + 1,nel
           do j2 = j1 + 1,nel
              if (.not.ormas_chkorth(i2,j2,istr,jstr,nela,orb)) cycle
              sgn12 = sgn1 * (-1)**(i2+j2)
              call util_window2(nel,i1,i2,j1,j2,s0,work)
              ovlp = ovlp + util_det(nel-2,thrdet,work)*sgn12
           end do
        end do
     end do
  end do

end subroutine ormas_str_dot_str_0_2
!################################################################################
subroutine ormas_str_dot_str_0_3(istr,jstr,ncore,nel,nela,orb,s0,work,ovlp)

  use,intrinsic :: iso_c_binding
  use mod_ormas,only : thrdet,nfun
  use mod_const,only : czero,runit

  implicit none
  !--------------------------------------------------------------------
  integer(c_int),intent(in) :: istr,jstr,ncore,nel,nela,orb(0:nela,*)
  complex(c_double_complex),intent(in) :: s0(1:nel,1:nel)
  complex(c_double_complex),intent(out) :: work(1:(nel-3),1:(nel-3))
  complex(c_double_complex),intent(out) :: ovlp
  !--------------------------------------------------------------------
  integer(c_int) :: i1,i2,i3,j1,j2,j3,sgn1,sgn12,sgn123
  complex(c_double_complex),external :: util_det
  logical(c_bool),external :: ormas_chkorth

  ovlp = czero
  do i1 = 1,nel
     do j1 = 1,nel
        if (.not.ormas_chkorth(i1,j1,istr,jstr,nela,orb)) cycle
        sgn1 = (-1)**(i1+j1)
        do i2 = i1 + 1,nel
           do j2 = j1 + 1,nel
              if (.not.ormas_chkorth(i2,j2,istr,jstr,nela,orb)) cycle
              sgn12 = sgn1 * (-1)**(i2+j2)
              do i3 = i2 + 1,nel
                 do j3 = j2 + 1,nel
                    if (.not.ormas_chkorth(i3,j3,istr,jstr,nela,orb)) cycle
                    sgn123 = sgn12 * (-1)**(i3+j3)
                    call util_window3(nel,i1,i2,i3,j1,j2,j3,s0,work)
                    ovlp = ovlp + util_det(nel-3,thrdet,work)*sgn123
                 end do
              end do
           end do
        end do
     end do
  end do

end subroutine ormas_str_dot_str_0_3
!################################################################################
subroutine ormas_str_dot_str_0_4(istr,jstr,ncore,nel,nela,orb,s0,work,ovlp)

  use,intrinsic :: iso_c_binding
  use mod_ormas,only : thrdet,nfun
  use mod_const,only : czero,runit

  implicit none
  !--------------------------------------------------------------------
  integer(c_int),intent(in) :: istr,jstr,ncore,nel,nela,orb(0:nela,*)
  complex(c_double_complex),intent(in) :: s0(1:nel,1:nel)
  complex(c_double_complex),intent(out) :: work(1:(nel-5),1:(nel-5))
  complex(c_double_complex),intent(out) :: ovlp
  !--------------------------------------------------------------------
  integer(c_int) :: i1,i2,i3,i4,j1,j2,j3,j4,sgn1,sgn12,sgn123,sgn1234
  complex(c_double_complex),external :: util_det
  logical(c_bool),external :: ormas_chkorth

  ovlp = czero
  do i1 = 1,nel
     do j1 = 1,nel
        if (.not.ormas_chkorth(i1,j1,istr,jstr,nela,orb)) cycle
        sgn1 = (-1)**(i1+j1)
        do i2 = i1 + 1,nel
           do j2 = j1 + 1,nel
              if (.not.ormas_chkorth(i2,j2,istr,jstr,nela,orb)) cycle
              sgn12 = sgn1 * (-1)**(i2+j2)
              do i3 = i2 + 1,nel
                 do j3 = j2 + 1,nel
                    if (.not.ormas_chkorth(i3,j3,istr,jstr,nela,orb)) cycle
                    sgn123 = sgn12 * (-1)**(i3+j3)
                    do i4 = i3 + 1,nel
                       do j4 = j3 + 1,nel
                          if (.not.ormas_chkorth(i4,j4,istr,jstr,nela,orb)) cycle
                          sgn1234 = sgn123 * (-1)**(i4+j4)
                          call util_window4(nel,i1,i2,i3,i4,j1,j2,j3,j4,s0,work)
                          ovlp = ovlp + util_det(nel-4,thrdet,work)*sgn1234
                       end do
                    end do
                 end do
              end do
           end do
        end do
     end do
  end do

end subroutine ormas_str_dot_str_0_4
!################################################################################
