!######################################################################
subroutine bas_gen_pmat(pmat)
!
  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad, ndvr, xrad, radp, ecs_flag, irad_ecs, theta, wrad, cwrad, mapb
  use mod_const, only : zero, iunit
!
! atomic hamiltonian in banded storage
!
! example: KL = 2, KU = 1:
! <--------- input ---------->     <--------- output --------->
!  *    *    *    +    +    +       *    *    *   u14  u25  u36
!  *    *    +    +    +    +       *    *   u13  u24  u35  u46
!  *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
! a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
! a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *
! a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *
!
  implicit none
  complex(c_double_complex), intent(out) :: pmat(1:(2*ndvr+1), 1:(nrad-1))
  integer(c_int) :: ifun, irad, jrad, jll, jul, ij, jb1, ctdep1, ctdep2

  pmat(1:(2*ndvr+1), 1:(nrad-1)) = zero
  
  if (ecs_flag == 0) then
     do irad = 1, nrad - 1
        jll = max(1,        irad - ndvr)
        jul = min(nrad - 1, irad + ndvr)
        do jrad = jll, jul
           jb1 = ndvr + 1 + jrad - irad
           pmat(jb1, irad) = radp(jb1, irad)
        end do
     end do
  end if
  
  if (ecs_flag == 1) then
     ctdep1 = ndvr * 2 + 1
     ctdep2 = ndvr
     do irad = 1, nrad - 1
        jll = max(1,        irad - ndvr)
        jul = min(nrad - 1, irad + ndvr)
        do jrad = jll, jul
           jb1 = ndvr + 1 + jrad - irad

           if( irad < irad_ecs  .and. irad + ndvr < irad_ecs) then 
              ! not ECS area
!!!!!!!!!debug
              !write(6, "('not ECS:', f20.10)") 
              !write(*, "('irad =', I5)") irad
              !write(*, "('jb1 =', I5)") jb1
              !write(*, "('xrad =', f20.10)") xrad(irad)
              !write(*, "('radk =',f20.10)") radk(jb1, irad)
              !write(*, "('l =',I5)") l
              !write(*, '(I10)') 1000000000000
              pmat(jb1, irad) = radp(jb1, irad)

           else if( irad < irad_ecs .and. irad + ndvr >= irad_ecs ) then

              if(mapb(irad) == 0) then !first finite element
                 if(jb1 == ctdep1 - 1) then
                    ! for ecs, when left boundary ife = I-1 (ECS will start from ife = I.)
                    ! In my note, blue area
!!!!!!!!!!!!!!debug
                    !write(6, "('blue:', f20.10)") 
                    !write(*, "('irad =', I5)") irad
                    !write(*, "('jb1 =', I5)") jb1
                    !write(*, "('xrad =', f20.10)") xrad(irad)
                    !write(*, "('radk =',f20.10)") radk(jb1, irad)
                    !write(*, "('l =',I5)") l
                    !write(*, '(I10)') 1000000000000
                    pmat(jb1, irad) = radp(jb1, irad) * sqrt( wrad(irad_ecs) / cwrad(irad_ecs) )
                    ctdep1 = ctdep1 - 1
                 else 
                    ! not ECS area
!!!!!!!!!!debug
                    !write(6, "('not ECS:', f20.10)") 
                    !write(*, "('irad =', I5)") irad
                    !write(*, "('jb1 =', I5)") jb1
                    !write(*, "('xrad =', f20.10)") xrad(irad)
                    !write(*, "('radk =',f20.10)") radk(jb1, irad)
                    !write(*, "('l =',I5)") l
                    !write(*, '(I10)') 1000000000000
                    pmat(jb1, irad) = radp(jb1, irad)
                 end if
              else !not first finite element
                 if(jb1 == ctdep1) then
                    ! for ecs, when left boundary ife = I-1 (ECS will start from ife = I.)
                    ! In my note, blue area
!!!!!!!!!!!debug
                    !write(6, "('blue:', f20.10)")  
                    !write(*, "('irad =', I5)") irad
                    !write(*, "('jb1 =', I5)") jb1
                    !write(*, "('xrad =', f20.10)") xrad(irad)
                    !write(*, "('radk =',f20.10)") radk(jb1, irad)
                    !write(*, "('l =',I5)") l
                    !write(*, "('cwrad(mapf(mapb(irad))+1) =',f20.10)") cwrad(mapf(mapb(irad)+1))
                    !write(*, '(I10)') 1000000000000
                    pmat(jb1, irad) = radp(jb1, irad) * sqrt( wrad(irad_ecs) / cwrad(irad_ecs) )
                    ctdep1 = ctdep1 - 1
                 else 
                    ! not ECS area
                    !write(6, "('not ECS:', f20.10)") 
                    !write(*, "('irad =', I5)") irad
                    !write(*, "('jb1 =', I5)") jb1
                    !write(*, "('xrad =', f20.10)") xrad(irad)
                    !write(*, "('radk =',f20.10)") radk(jb1, irad)
                    !write(*, "('l =',I5)") l
                    !write(*, '(I10)') 1000000000000
                    pmat(jb1, irad) = radp(jb1, irad)
                 end if
              end if

           else if( irad == irad_ecs ) then
              if(jb1 < ndvr + 1) then
                 ! for ecs, when upper boundary between ife = I - 1 and I  (ECS will start from ife = I.)
                 ! In my note, blue area
!!!!!!!!!!debug
                 !write(6, "('blue:', f20.10)") 
                 !write(*, "('irad =', I5)") irad
                 !write(*, "('jb1 =', I5)") jb1
                 !write(*, "('xrad =', f20.10)") xrad(irad)
                 !write(*, "('radk =',f20.10)") radk(jb1, irad)
                 !write(*, "('l =',I5)") l
                 !write(*, '(I10)') 1000000000000
                 pmat(jb1, irad) = radp(jb1, irad) * sqrt( wrad(irad) / cwrad(irad) )
              else if(jb1 == ndvr + 1) then
                 ! for ecs, when center of boundary between ife = I - 1 and I  (ECS will start from ife = I.)
                 ! In my note, black area
!!!!!!!!!!!!!!!!debug
                 !write(6, "('++++++++black+++++++:', f20.10)") 
                 !write(*, "('irad =', I5)") irad
                 !write(*, "('jb1 =', I5)") jb1
                 !write(*, "('xrad =', f20.10)") xrad(irad)
                 !write(*, "('radk =',f20.10)") radk(jb1, irad)
                 !write(*, "('l =',I5)") l
                 !write(*, '(I10)') 1000000000000
                 pmat(jb1, irad) = radp(jb1, irad) * wrad(irad) / cwrad(irad)
              else if(jb1 > ndvr + 1) then
                 ! for ecs, when lower boundary between ife = I - 1 and I  (ECS will start from ife = I.)
                 ! In my note, red area
!!!!!!!!!!!!!!!!debug
                 !write(6, "('red:', f20.10)") 
                 !write(*, "('irad =', I5)") irad
                 !write(*, "('jb1 =', I5)") jb1
                 !write(*, "('xrad =', f20.10)") xrad(irad)
                 !write(*, "('radk =',f20.10)") radk(jb1, irad)
                 !write(*, "('l =',I5)") l
                 !write(*, '(I10)') 1000000000000
                 pmat(jb1, irad) = radp(jb1, irad) * sqrt( wrad(irad) / cwrad(irad) ) * EXP((-0.5) * iunit * theta)
              end if
              
           else if( irad > irad_ecs .and. irad - ndvr <= irad_ecs ) then
              if(jb1 == ctdep2) then
                 ! for ecs, when right boundary between ife = I  (ECS will start from ife = I.)
                 ! In my note, red area
!!!!!!!!!!!!!!!!debug
                 !write(6, "('red:', f20.10)")  
                 !write(*, "('irad =', I5)") irad
                 !write(*, "('jb1 =', I5)") jb1
                 !write(*, "('xrad =', f20.10)") xrad(irad)
                 !write(*, "('radk =',f20.10)") radk(jb1, irad)
                 !write(*, "('cwrad(mapf(mapb(irad))) =',f20.10)") cwrad(irad_ecs)
                 !write(*, "('l =',I5)") l
                 !write(*, '(I10)') 1000000000000
                 pmat(jb1, irad) = radp(jb1, irad) * sqrt( wrad(irad_ecs) / cwrad(irad_ecs) ) * &
                      & EXP( (-0.5) * iunit * theta)
                 ctdep2 = ctdep2 - 1
              else
!!!!!!!!!!!!!!!!debug
                 !write(6, "('C.ECS:', f20.10)")  
                 !write(*, "('irad =', I5)") irad
                 !write(*, "('jb1 =', I5)") jb1
                 !!write(*, "('xrad =', f20.10)") xrad(irad)        
                 !write(*, "('radk =',f20.10)") radk(jb1, irad)
                 !write(*, "('l =',I5)") l
                 !write(*, '(I10)') 1000000000000
                 pmat(jb1, irad) = radp(jb1, irad) * EXP( (-1.0) * iunit * theta)
              end if

           else if( irad - ndvr > irad_ecs  .and. irad - ndvr >  irad_ecs) then
              ! completely ECS area
!!!!!!!!!!!!!!!!debug
              !write(6, "('C.ECS:', f20.10)")  
              !write(*, "('irad =', I5)") irad
              !write(*, "('jb1 =', I5)") jb1
              !write(*, "('xrad =', f20.10)") xrad(irad)
              !write(*, "('radk =',f20.10)") radk(jb1, irad)
              !write(*, "('l =',I5)") l
              !write(*, '(I10)') 1000000000000
              pmat(jb1, irad) = radp(jb1, irad) * EXP((-1.0) * iunit * theta)   
           end if
           
        end do
     end do
  end if


end subroutine bas_gen_pmat
!######################################################################
