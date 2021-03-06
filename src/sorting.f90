!===============================================================================!
MODULE MOD_Sorting
!-------------------------------------------------------------------------------!
! Sort a set of numbers.                                                        !
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
PRIVATE
!-------------------------------------------------------------------------------!
INTERFACE SortArray
  MODULE PROCEDURE SortArray
END INTERFACE
!-------------------------------------------------------------------------------!
PUBLIC :: SortArray
!-------------------------------------------------------------------------------!
!
!
!
!===============================================================================!
CONTAINS
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE SortArray(x)
!-------------------------------------------------------------------------------!
! Receives an array x and sorts it into ascending order                         !
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(INOUT) :: x(1:)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER            :: i
INTEGER            :: Location
!-------------------------------------------------------------------------------!

DO i = 1, SIZE(x)-1                     ! except for the last
  Location = FindMinimum(x,i,SIZE(x))   ! find min from this to last
  CALL Swap(x(i),x(Location))           ! swap this and the minimum
END DO
   
!-------------------------------------------------------------------------------!
END SUBROUTINE SortArray
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE Swap(a,b)
!-------------------------------------------------------------------------------!
! Swaps the values of its two formal arguments                                  !
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(INOUT) :: a, b
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL               :: temp
!-------------------------------------------------------------------------------!

temp = a
a    = b
b    = temp

!-------------------------------------------------------------------------------!
END SUBROUTINE Swap
!===============================================================================!
!
!
!
!===============================================================================!
FUNCTION FindMinimum(x,StartPosition,EndPosition)
!-------------------------------------------------------------------------------!
! Returns the location of the minimum in the section between Start and End      !
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,   INTENT(IN) :: x(1:)
INTEGER,INTENT(IN) :: StartPosition, EndPosition
INTEGER            :: FindMinimum
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL               :: Minimum
INTEGER            :: Location
INTEGER            :: i
!-------------------------------------------------------------------------------!

Minimum  = x(StartPosition)           ! assume the first is the min
Location = StartPosition              ! record its position
DO i = StartPosition+1, EndPosition   ! start with next elements
    IF (x(i) < Minimum) THEN          ! if x(i) less than the min?
      Minimum  = x(i)                 ! Yes, a new minimum found
      Location = i                    ! record its position
    END IF
END DO
FindMinimum = Location                ! return the position

!-------------------------------------------------------------------------------!
END FUNCTION  FindMinimum
!===============================================================================!
!
!
!
!===============================================================================!
END MODULE MOD_Sorting
!===============================================================================!
