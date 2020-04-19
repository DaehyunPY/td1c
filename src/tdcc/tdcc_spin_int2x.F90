!###################################################################################
integer function tdcc_spin_int2x(s1, s2, s3, s4)

!  1 : <p1 q1||r1 s1> = +(pq|rs)-(pq|sr)
!  2 : <p2 q2||r2 s2> = +(pq|rs)-(pq|sr) = 1
!  3 : <p1 q2||r1 s2> = +(pq|rs)
!  4 : <p2 q1||r2 s1> = +(pq|rs)         = 3
!  5 : <p1 q2||r2 s1> = -(pq|sr)
!  6 : <p2 q1||r1 s2> = -(pq|sr)         = 6
!  0 : others

  implicit none
  integer, intent(in) :: s1, s2, s3, s4
  integer :: spin_type

  if (s1 + s2 .ne. s3 + s4) then
     spin_type = 0
  else if (s1 == 1 .and. s2 == 1 .and. s3 == 1 .and. s4 == 1) then
     spin_type = 1
  else if (s1 == 2 .and. s2 == 2 .and. s3 == 2 .and. s4 == 2) then
     spin_type = 2
  else if (s1 == 1 .and. s2 == 2 .and. s3 == 1 .and. s4 == 2) then
     spin_type = 3
  else if (s1 == 2 .and. s2 == 1 .and. s3 == 2 .and. s4 == 1) then
     spin_type = 4
  else if (s1 == 1 .and. s2 == 2 .and. s3 == 2 .and. s4 == 1) then
     spin_type = 5
  else if (s1 == 2 .and. s2 == 1 .and. s3 == 1 .and. s4 == 2) then
     spin_type = 6
  else
     spin_type = 0
  end if

  tdcc_spin_int2x = spin_type

end function tdcc_spin_int2x
!###################################################################################
