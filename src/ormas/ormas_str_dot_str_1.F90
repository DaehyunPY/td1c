!################################################################################
subroutine ormas_str_dot_str_1(iipx,istr,jstr,ncore,nel,nela,orb,s0,work,dmat)

  use,intrinsic :: iso_c_binding
  use mod_ormas,only : thrdet,nfun,nfcore
  use mod_const,only : czero,runit

  implicit none
  !--------------------------------------------------------------------
  integer(c_int),intent(in) :: iipx,istr,jstr,ncore,nel,nela,orb(0:nela,*)
  complex(c_double_complex),intent(in) :: s0(1:nel,1:nel)
  complex(c_double_complex),intent(out) :: work(1:(nel-iipx),1:(nel-iipx))
  complex(c_double_complex),intent(out) :: dmat(1:nfun,1:nfun)
  !--------------------------------------------------------------------

  if (iipx == 0) then
     call ormas_str_dot_str_1_0(istr,jstr,ncore,nel,nela,orb,s0,work,dmat)
  else if (iipx == 1) then
     call ormas_str_dot_str_1_1(istr,jstr,ncore,nel,nela,orb,s0,work,dmat)
  else if (iipx == 2) then
     call ormas_str_dot_str_1_2(istr,jstr,ncore,nel,nela,orb,s0,work,dmat)
  else if (iipx == 3) then
     call ormas_str_dot_str_1_3(istr,jstr,ncore,nel,nela,orb,s0,work,dmat)
  else if (iipx == 4) then
     call ormas_str_dot_str_1_4(istr,jstr,ncore,nel,nela,orb,s0,work,dmat)
  else
     stop 'ormas_str_dot_str_1: bad iipx.'
  end if

  return

end subroutine ormas_str_dot_str_1
!################################################################################
subroutine ormas_str_dot_str_1_0(istr,jstr,ncore,nel,nela,orb,s0,work,dmat)

  use,intrinsic :: iso_c_binding
  use mod_ormas,only : thrdet,nfun,nfcore
  use mod_const,only : czero,runit

  implicit none
  !--------------------------------------------------------------------
  integer(c_int),intent(in) :: istr,jstr,ncore,nel,nela,orb(0:nela,*)
  complex(c_double_complex),intent(in) :: s0(1:nel,1:nel)
  complex(c_double_complex),intent(out) :: work(1:(nel-1),1:(nel-1))
  complex(c_double_complex),intent(out) :: dmat(1:nfun,1:nfun)
  !--------------------------------------------------------------------
  integer(c_int) :: ifun,jfun,i0,j0,sgn0
  complex(c_double_complex),external :: util_det
  integer(c_int),external :: ormas_elec2fun

  dmat = czero
  do i0 = 1,nel
     ifun = ormas_elec2fun(i0,istr,nela,orb)
     do j0 = 1,nel
        jfun = ormas_elec2fun(j0,jstr,nela,orb)
        sgn0 = (-1)**(i0+j0)
        call util_window1(nel,i0,j0,s0,work)
        dmat(ifun,jfun) = dmat(ifun,jfun) + util_det(nel-1,thrdet,work) * sgn0
     end do
  end do

end subroutine ormas_str_dot_str_1_0
!################################################################################
subroutine ormas_str_dot_str_1_1(istr,jstr,ncore,nel,nela,orb,s0,work,dmat)

  use,intrinsic :: iso_c_binding
  use mod_ormas,only : thrdet,nfun,nfcore
  use mod_const,only : czero,runit

  implicit none
  !--------------------------------------------------------------------
  integer(c_int),intent(in) :: istr,jstr,ncore,nel,nela,orb(0:nela,*)
  complex(c_double_complex),intent(in) :: s0(1:nel,1:nel)
  complex(c_double_complex),intent(out) :: work(1:(nel-2),1:(nel-2))
  complex(c_double_complex),intent(out) :: dmat(1:nfun,1:nfun)
  !--------------------------------------------------------------------
  integer(c_int) :: ifun,jfun,i0,i1,j0,j1,i0x,i1x,j0x,j1x,sgn0,sgni,sgnj,sgn01
  complex(c_double_complex),external :: util_det
  integer(c_int),external :: ormas_elec2fun
  logical(c_bool),external :: ormas_chkorth

  dmat = czero
  do i0 = 1,nel
     ifun = ormas_elec2fun(i0,istr,nela,orb)
     do j0 = 1,nel
        jfun = ormas_elec2fun(j0,jstr,nela,orb)
        sgn0 = (-1)**(i0+j0)
        do i1 = 1,nel
           if (i1 == i0) cycle
           call ormas_sgnperm2(i0,i1,i0x,i1x,sgni)
           do j1 = 1,nel
              if (j1 == j0) cycle
              if (.not.ormas_chkorth(i1,j1,istr,jstr,nela,orb)) cycle
              call ormas_sgnperm2(j0,j1,j0x,j1x,sgnj)
              sgn01 = sgn0*sgni*sgnj
              call util_window2(nel,i0x,i1x,j0x,j1x,s0,work)
              dmat(ifun,jfun) = dmat(ifun,jfun) + util_det(nel-2,thrdet,work)*sgn01
           end do
        end do
     end do
  end do

end subroutine ormas_str_dot_str_1_1
!################################################################################
subroutine ormas_str_dot_str_1_2(istr,jstr,ncore,nel,nela,orb,s0,work,dmat)

  use,intrinsic :: iso_c_binding
  use mod_ormas,only : thrdet,nfun,nfcore
  use mod_const,only : czero,runit

  implicit none
  !--------------------------------------------------------------------
  integer(c_int),intent(in) :: istr,jstr,ncore,nel,nela,orb(0:nela,*)
  complex(c_double_complex),intent(in) :: s0(1:nel,1:nel)
  complex(c_double_complex),intent(out) :: work(1:(nel-3),1:(nel-3))
  complex(c_double_complex),intent(out) :: dmat(1:nfun,1:nfun)
  !--------------------------------------------------------------------
  integer(c_int) :: ifun,jfun,i0,i1,i2,j0,j1,j2,i0x,i1x,i2x,j0x,j1x,j2x,&
       sgn0,sgni,sgnj,sgn012
  complex(c_double_complex),external :: util_det
  integer(c_int),external :: ormas_elec2fun
  logical(c_bool),external :: ormas_chkorth

  dmat = czero
  do i0 = 1,nel
     ifun = ormas_elec2fun(i0,istr,nela,orb)
     do j0 = 1,nel
        jfun = ormas_elec2fun(j0,jstr,nela,orb)
        sgn0 = (-1)**(i0+j0)
        do i1 = 1,nel
           if (i1 == i0) cycle
           do i2 = i1 + 1,nel
              if (i2 == i0) cycle
              call ormas_sgnperm3(i0,i1,i2,i0x,i1x,i2x,sgni)
              do j1 = 1,nel
                 if (j1 == j0) cycle
                 if (.not.ormas_chkorth(i1,j1,istr,jstr,nela,orb)) cycle
                 do j2 = j1 + 1,nel
                    if (j2 == j0) cycle
                    if (.not.ormas_chkorth(i2,j2,istr,jstr,nela,orb)) cycle
                    call ormas_sgnperm3(j0,j1,j2,j0x,j1x,j2x,sgnj)
                    sgn012 = sgn0*sgni*sgnj
                    call util_window3(nel,i0x,i1x,i2x,j0x,j1x,j2x,s0,work)
                    dmat(ifun,jfun) = dmat(ifun,jfun) + util_det(nel-3,thrdet,work)*sgn012
                 end do
              end do
           end do
        end do
     end do
  end do

end subroutine ormas_str_dot_str_1_2
!################################################################################
subroutine ormas_str_dot_str_1_3(istr,jstr,ncore,nel,nela,orb,s0,work,dmat)

  use,intrinsic :: iso_c_binding
  use mod_ormas,only : thrdet,nfun,nfcore
  use mod_const,only : czero,runit

  implicit none
  !--------------------------------------------------------------------
  integer(c_int),intent(in) :: istr,jstr,ncore,nel,nela,orb(0:nela,*)
  complex(c_double_complex),intent(in) :: s0(1:nel,1:nel)
  complex(c_double_complex),intent(out) :: work(1:(nel-4),1:(nel-4))
  complex(c_double_complex),intent(out) :: dmat(1:nfun,1:nfun)
  !--------------------------------------------------------------------
  integer(c_int) :: ifun,jfun,i0,i1,i2,i3,j0,j1,j2,j3,i0x,i1x,i2x,i3x,&
       j0x,j1x,j2x,j3x,sgn0,sgni,sgnj,sgn0123
  complex(c_double_complex),external :: util_det
  integer(c_int),external :: ormas_elec2fun
  logical(c_bool),external :: ormas_chkorth

  dmat = czero
  do i0 = 1,nel
     ifun = ormas_elec2fun(i0,istr,nela,orb)
     do j0 = 1,nel
        jfun = ormas_elec2fun(j0,jstr,nela,orb)
        sgn0 = (-1)**(i0+j0)
        do i1 = 1,nel
           if (i1 == i0) cycle
           do i2 = i1 + 1,nel
              if (i2 == i0) cycle
              do i3 = i2 + 1,nel
                 if (i3 == i0) cycle
                 call ormas_sgnperm4(i0,i1,i2,i3,i0x,i1x,i2x,i3x,sgni)
                 do j1 = 1,nel
                    if (j1 == j0) cycle
                    if (.not.ormas_chkorth(i1,j1,istr,jstr,nela,orb)) cycle
                    do j2 = j1 + 1,nel
                       if (j2 == j0) cycle
                       if (.not.ormas_chkorth(i2,j2,istr,jstr,nela,orb)) cycle
                       do j3 = j2 + 1,nel
                          if (j3 == j0) cycle
                          if (.not.ormas_chkorth(i3,j3,istr,jstr,nela,orb)) cycle
                          call ormas_sgnperm4(j0,j1,j2,j3,j0x,j1x,j2x,j3x,sgnj)
                          sgn0123 = sgn0*sgni*sgnj
                          call util_window4(nel,i0x,i1x,i2x,i3x,j0x,j1x,j2x,j3x,s0,work)
                          dmat(ifun,jfun) = dmat(ifun,jfun) + util_det(nel-4,thrdet,work)*sgn0123
                       end do
                    end do
                 end do
              end do
           end do
        end do
     end do
  end do

end subroutine ormas_str_dot_str_1_3
!################################################################################
subroutine ormas_str_dot_str_1_4(istr,jstr,ncore,nel,nela,orb,s0,work,dmat)

  use,intrinsic :: iso_c_binding
  use mod_ormas,only : thrdet,nfun,nfcore
  use mod_const,only : czero,runit

  implicit none
  !--------------------------------------------------------------------
  integer(c_int),intent(in) :: istr,jstr,ncore,nel,nela,orb(0:nela,*)
  complex(c_double_complex),intent(in) :: s0(1:nel,1:nel)
  complex(c_double_complex),intent(out) :: work(1:(nel-5),1:(nel-5))
  complex(c_double_complex),intent(out) :: dmat(1:nfun,1:nfun)
  !--------------------------------------------------------------------
  integer(c_int) :: ifun,jfun,i0,i1,i2,i3,i4,j0,j1,j2,j3,j4,i0x,i1x,i2x,i3x,i4x,&
       j0x,j1x,j2x,j3x,j4x,sgn0,sgni,sgnj,sgn01234
  complex(c_double_complex),external :: util_det
  integer(c_int),external :: ormas_elec2fun
  logical(c_bool),external :: ormas_chkorth

  dmat = czero
  do i0 = 1,nel
     ifun = ormas_elec2fun(i0,istr,nela,orb)
     do j0 = 1,nel
        jfun = ormas_elec2fun(j0,jstr,nela,orb)
        sgn0 = (-1)**(i0+j0)
        do i1 = 1,nel
           if (i1 == i0) cycle
           do i2 = i1 + 1,nel
              if (i2 == i0) cycle
              do i3 = i2 + 1,nel
                 if (i3 == i0) cycle
                 do i4 = i3 + 1,nel
                    if (i4 == i0) cycle
                    call ormas_sgnperm5(i0,i1,i2,i3,i4,i0x,i1x,i2x,i3x,i4x,sgni)
                    do j1 = 1,nel
                       if (j1 == j0) cycle
                       if (.not.ormas_chkorth(i1,j1,istr,jstr,nela,orb)) cycle
                       do j2 = j1 + 1,nel
                          if (j2 == j0) cycle
                          if (.not.ormas_chkorth(i2,j2,istr,jstr,nela,orb)) cycle
                          do j3 = j2 + 1,nel
                             if (j3 == j0) cycle
                             if (.not.ormas_chkorth(i3,j3,istr,jstr,nela,orb)) cycle
                             do j4 = j3 + 1,nel
                                if (j4 == j0) cycle
                                if (.not.ormas_chkorth(i4,j4,istr,jstr,nela,orb)) cycle
                                call ormas_sgnperm5(j0,j1,j2,j3,j4,j0x,j1x,j2x,j3x,j4x,sgnj)
                                sgn01234 = sgn0*sgni*sgnj
                                call util_window5(nel,i0x,i1x,i2x,i3x,i4x,j0x,j1x,j2x,j3x,j4x,s0,work)
                                dmat(ifun,jfun) = dmat(ifun,jfun) + util_det(nel-5,thrdet,work)*sgn01234
                             end do
                          end do
                       end do
                    end do
                 end do
              end do
           end do
        end do
     end do
  end do

end subroutine ormas_str_dot_str_1_4
!################################################################################
