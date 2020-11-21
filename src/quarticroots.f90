!===============================================================================!
MODULE MOD_QuarticRoots
!-------------------------------------------------------------------------------!
! Solve for the roots of a polynomial equation with real coefficients,          !
! up to quartic order.                                                          !
! Retrun a code indicating the nature of the roots found.                       !
!                                                                               !
! AUTHORS                                                                       !
!   - Alfred Morris,      Naval Surface Weapons Center, Dahlgren,VA             !
!   - William Davis,      Naval Surface Weapons Center, Dahlgren,VA             !
!   - Alan Miller,        CSIRO Mathematical & Information Sciences             !
!                         CLAYTON, VICTORIA, AUSTRALIA 3169                     !
!                         http://www.mel.dms.csiro.au/~alan                     !
!   - Ralph Carmichael,   Public Domain Aeronautical Software                   !
!                         http://www.pdas.com                                   !
! REVISION HISTORY                                                              !
!   DATE  VERS PERSON  STATEMENT OF CHANGES                                     !
!    ??    1.0 AHM&WLD Original coding                                          !
! 27Feb97  1.1   AM    Converted to be compatible with ELF90                    !
! 12Jul98  1.2   RLC   Module format; numerous style changes                    !
!  4Jan99  1.3   RLC   Made the tests for zero constant term exactly zero       !
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
PRIVATE
!-------------------------------------------------------------------------------!
INTEGER,       PARAMETER :: DP  = KIND(1.0D+00)
REAL(KIND=DP), PARAMETER :: EPS = EPSILON(1.0D+00)
!-------------------------------------------------------------------------------!
INTERFACE QuarticRoots
  MODULE PROCEDURE QuarticRoots
END INTERFACE
!-------------------------------------------------------------------------------!
PUBLIC :: QuarticRoots
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
FUNCTION CUBEROOT(x) RESULT(f)
!-------------------------------------------------------------------------------!
! Compute the Cube Root of a real number.                                       !
! If the argument is negative, then the cube root is also negative.             !
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL(KIND=DP),INTENT(IN) :: x
REAL(KIND=DP)            :: f
!-------------------------------------------------------------------------------!

IF (x .LT. 0.0) THEN
  f = -EXP(LOG(-x)/3.0)
ELSE IF (x .GT. 0.0) THEN
  f = EXP(LOG(x)/3.0)
ELSE
  f = 0.0
END IF
RETURN

!-------------------------------------------------------------------------------!
END FUNCTION CUBEROOT
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE LinearRoot(a, z)
!-------------------------------------------------------------------------------!
! Computes the roots of the real polynomial                                     !
!   a(1) + a(2)*z                                                               !
! and stores the results in z. It is assumed that a(2) is non-zero.             !
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL(KIND=DP),INTENT(IN)  :: a(:)
REAL(KIND=DP),INTENT(OUT) :: z
!-------------------------------------------------------------------------------!

IF (a(2) .EQ. 0.0) THEN
  z = 0.0
ELSE
  z = -a(1)/a(2)
END IF
RETURN

!-------------------------------------------------------------------------------!
END SUBROUTINE LinearRoot
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE OneLargeTwoSmall(a1, a2, a4, w, z)
!-------------------------------------------------------------------------------!
! Compute the roots of a cubic when one root,                                   !
! w, is known to be much larger in magnitude than the other two                 !
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL(KIND=DP),   INTENT(IN)  :: a1, a2, a4
REAL(KIND=DP),   INTENT(IN)  :: w
COMPLEX(KIND=DP),INTENT(OUT) :: z(:)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL(KIND=DP)                :: aq(3)
!-------------------------------------------------------------------------------!

aq(1) = a1
aq(2) = a2+a1/w
aq(3) = -a4*w

CALL QuadraticRoots(aq, z)

z(3) = CMPLX(w,0.0,DP)

IF (AIMAG(z(1)) .EQ. 0.0) RETURN

z(3) = z(2)
z(2) = z(1)
z(1) = CMPLX(w,0.0,DP)
RETURN

!-------------------------------------------------------------------------------!
END SUBROUTINE OneLargeTwoSmall
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE QuadraticRoots(a, z)
!-------------------------------------------------------------------------------!
! Computes the roots of the real polynomial                                     !
!   a(1) + a(2)*z + a(3)*z**2                                                   !
! and stores the results in z. It is assumed that a(3) is non-zero.             !
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL(KIND=DP),   INTENT(IN)  :: a(:)
COMPLEX(KIND=DP),INTENT(OUT) :: z(:)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL(KIND=DP)                :: d, r, w, x, y
!-------------------------------------------------------------------------------!

IF (a(1) .EQ. 0.0) THEN           ! EPS is a global module constant (private)
  z(1) = (0.0, 0.0)               ! one root is obviously 0.0
  z(2) = CMPLX(-a(2)/a(3), 0.0,DP)! remainder is a linear eq.
  RETURN
END IF

d = a(2)*a(2) - 4.0*a(1)*a(3)     ! the discriminant
IF (ABS(d) .LE. 2.0*EPS*a(2)*a(2)) THEN
  z(1) = CMPLX(-0.5*a(2)/a(3), 0.0, DP)   ! discriminant is tiny
  z(2) = z(1)
  RETURN
END IF

r = SQRT(ABS(d))
IF (d .LT. 0.0) THEN
  x    = -0.5*a(2)/a(3)           ! negative discriminant => roots are complex
  y    = ABS(0.5*r/a(3))
  z(1) = CMPLX(x, y, DP)
  z(2) = CMPLX(x,-y, DP)          ! its conjugate
  RETURN
END IF

IF (a(2) .NE. 0.0) THEN           ! see Numerical Recipes, sec. 5.5
  w    = -(a(2) + SIGN(r,a(2)))
  z(1) = CMPLX(2.0*a(1)/w,  0.0, DP)
  z(2) = CMPLX(0.5*w/a(3), 0.0, DP)
  RETURN
END IF

x = ABS(0.5*r/a(3))               ! a(2)=0 if you get here
z(1) = CMPLX( x, 0.0, DP)
z(2) = CMPLX(-x, 0.0, DP)
RETURN

!-------------------------------------------------------------------------------!
END SUBROUTINE QuadraticRoots
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE CubicRoots(a, z)
!-------------------------------------------------------------------------------!
! Compute the roots of the real polynomial                                      !
!   a(1) + a(2)*z + a(3)*z**2 + a(4)*z**3                                       !
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL(KIND=DP),   INTENT(IN)  :: a(:)
COMPLEX(KIND=DP),INTENT(OUT) :: z(:)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL(KIND=DP)                :: aq(3), arg, c, cf, d, p, p1, q, q1
REAL(KIND=DP)                :: r, ra, rb, rq, rt
REAL(KIND=DP)                :: r1, s, sf, sq, sum, t, tol, t1, w
REAL(KIND=DP)                :: w1, w2, x, x1, x2, x3, y, y1, y2, y3
!-------------------------------------------------------------------------------!

! It is assumed that a(4) is non-zero. No test is made here.

IF (a(1) .EQ. 0.0) THEN
  z(1) = (0.0, 0.0)                     ! one root is obviously zero
  CALL QuadraticRoots(a(2:4), z(2:3))   ! remaining 2 roots here
  RETURN
END IF

p   = a(3)/(3.0*a(4))
q   = a(2)/a(4)
r   = a(1)/a(4)
tol = 4.0*EPS

c = 0.0
t = a(2) - p*a(3)
IF (ABS(t) .GT. tol*ABS(a(2))) c = t/a(4)

t = 2.0*p*p - q
IF (ABS(t) .LE. tol*ABS(q)) t = 0.0

d = r + p*t
IF (ABS(d) .LE. tol*ABS(r)) GO TO 110

! Set  sq = (a(4)/s)**2 * (c**3/27 + d**2/4)

s  = MAX(ABS(a(1)), ABS(a(2)), ABS(a(3)))
p1 = a(3)/(3.0*s)
q1 = a(2)/s
r1 = a(1)/s

t1 = q - 2.25*p*p
IF (ABS(t1) .LE. tol*ABS(q)) t1 = 0.0

w  = 0.25*r1*r1
w1 = 0.5*p1*r1*t
w2 = q1*q1*t1/27.0

IF (w1 .GE. 0.0) THEN
  w  = w + w1
  sq = w + w2
ELSE IF (w2 .LT. 0.0) THEN
  sq = w + (w1 + w2)
ELSE
  w  = w + w2
  sq = w + w1
END IF

IF (ABS(sq) .LE. tol*w) sq = 0.0

rq = ABS(s/a(4))*SQRT(ABS(sq))
IF (sq .GE. 0.0) GO TO 40

! All roots are real

arg = ATAN2(rq, -0.5*d)
cf  = COS(arg/3.0)
sf  = SIN(arg/3.0)
rt  = SQRT(-c/3.0)
y1  = 2.0*rt*cf
y2  = -rt*(cf + SQRT(3.0)*sf)
y3  = -(d/y1)/y2

x1 = y1 - p
x2 = y2 - p
x3 = y3 - p

IF (ABS(x1) .GT. ABS(x2)) CALL SWAP(x1,x2)
IF (ABS(x2) .GT. ABS(x3)) CALL SWAP(x2,x3)
IF (ABS(x1) .GT. ABS(x2)) CALL SWAP(x1,x2)

w = x3

IF (ABS(x2) .LT. 0.1*ABS(x3)) GO TO 70
IF (ABS(x1) .LT. 0.1*ABS(x2)) x1 = - (r/x3)/x2
z(1) = CMPLX(x1, 0.0,DP)
z(2) = CMPLX(x2, 0.0,DP)
z(3) = CMPLX(x3, 0.0,DP)
RETURN

! Real and complex roots

40 ra = CUBEROOT(-0.5*d - SIGN(rq,d))
rb = -c/(3.0*ra)
t  = ra + rb
w  = -p
x  = -p

IF (ABS(t) .LE. tol*ABS(ra)) GO TO 41

w = t - p
x = -0.5*t - p
IF (ABS(x) .LE. tol*ABS(p)) x = 0.0

41 t = ABS(ra - rb)
y = 0.5*SQRT(3.0)*t
IF (t .LE. tol*ABS(ra)) GO TO 60
IF (ABS(x) .LT. ABS(y)) GO TO 50

s = ABS(x)
t = y/x
GO TO 51

50 s = ABS(y)
t = x/y
51 IF (s .LT. 0.1*ABS(w)) GO TO 70

w1 = w/s
sum = 1.0 + t*t
IF (w1*w1 .LT. 0.01*sum) w = - ((r/sum)/s)/s

z(1) = CMPLX(w, 0.0,DP)
z(2) = CMPLX(x, y,DP)
z(3) = CMPLX(x,-y,DP)
RETURN

! At least two roots are equal

60 IF (ABS(x) .LT. ABS(w)) GO TO 61
IF (ABS(w) .LT. 0.1*ABS(x)) w = - (r/x)/x

z(1) = CMPLX(w, 0.0,DP)
z(2) = CMPLX(x, 0.0,DP)
z(3) = z(2)
RETURN

61 IF (ABS(x) .LT. 0.1*ABS(w)) GO TO 70

z(1) = CMPLX(x, 0.0,DP)
z(2) = z(1)
z(3) = CMPLX(w, 0.0,DP)
RETURN

!-------------------------------------------------------------!
! Here w is much larger in magnitude than the other roots.    !
! as a result, the other roots may be exceedingly inaccurate  !
! because of roundoff error. To deal with this, a quadratic   !
! is formed whose roots are the same as the smaller roots of  !
! the cubic.  this quadratic is then solved.                  !
!-------------------------------------------------------------!
! This code was written by william l. davis (nswc).           !
!-------------------------------------------------------------!

70 aq(1) = a(1)
aq(2) = a(2) + a(1)/w
aq(3) = -a(4)*w
CALL QuadraticRoots(aq, z)
z(3) = CMPLX(w, 0.0,DP)

IF (AIMAG(z(1)) == 0.0) RETURN

z(3) = z(2)
z(2) = z(1)
z(1) = CMPLX(w, 0.0,DP)
RETURN
!-----------------------------------------------------------------------


!                   CASE WHEN D = 0

110 z(1) = CMPLX(-p, 0.0,DP)
  w = SQRT(ABS(c))
  IF (c .LT. 0.0) GO TO 120
  z(2) = CMPLX(-p, w,DP)
  z(3) = CMPLX(-p,-w,DP)
  RETURN

120 IF (p /= 0.0) GO TO 130
  z(2) = CMPLX(w, 0.0,DP)
  z(3) = CMPLX(-w, 0.0,DP)
  RETURN

130 x = -(p + SIGN(w,p))
  z(3) = CMPLX(x, 0.0,DP)
  t = 3.0*a(1)/(a(3)*x)
  IF (ABS(p) .GT. ABS(t)) GO TO 131
  z(2) = CMPLX(t, 0.0,DP)
  RETURN
131 z(2) = z(1)
  z(1) = CMPLX(t, 0.0,DP)
  RETURN

!-------------------------------------------------------------------------------!
END SUBROUTINE CubicRoots
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE QuarticRoots(a,z)
!-------------------------------------------------------------------------------!
! Compute the roots of the real polynomial                                      !
!   a(1) + a(2)*z + ... + a(5)*z**4                                             !
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL(KIND=DP),   INTENT(IN)  :: a(:)
COMPLEX(KIND=DP),INTENT(OUT) :: z(:)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
COMPLEX(KIND=DP)             :: w
REAL(KIND=DP)                :: b,b2, c, d, e, h, p, q, r, t
REAL(KIND=DP)                :: temp(4)
REAL(KIND=DP)                :: u, v, v1, v2, x, x1, x2, x3, y
!-------------------------------------------------------------------------------!

! NOTE - It is assumed that a(5) is non-0.0. No test is made here

IF (a(1) .EQ. 0.0) THEN
  z(1) = (0.0, 0.0)    !  1.0 root is obviously 0.0
  CALL CubicRoots(a(2:), z(2:))
  RETURN
END IF

b  = a(4)/(4.0*a(5))
c  = a(3)/a(5)
d  = a(2)/a(5)
e  = a(1)/a(5)
b2 = b*b

p = 0.5*(c - 6.0*b2)
q = d - 2.0*b*(c - 4.0*b2)
r = b2*(c - 3.0*b2) - b*d + e

!-----------------------------------------------------------------!
! Solve the resolvent cubic equation. The cubic has at least one  !
! nonnegative real root. If w1, w2, w3 are the roots of the cubic !
! then the roots of the originial equation are                    !
!   z = -b + csqrt(w1) + csqrt(w2) + csqrt(w3)                    !
! where the signs of the square roots are chosen so               !
! that csqrt(w1) * csqrt(w2) * csqrt(w3) = -q/8.                  !
!-----------------------------------------------------------------!

temp(1) = -q*q/64.0
temp(2) = 0.25*(p*p - r)
temp(3) =  p
temp(4) = 1.0
CALL CubicRoots(temp,z)
IF (AIMAG(z(2)) /= 0.0) GO TO 60

! The resolvent cubic has only real roots
! reorder the roots in increasing order

x1 = DBLE(z(1))
x2 = DBLE(z(2))
x3 = DBLE(z(3))
IF (x1 .GT. x2) CALL SWAP(x1,x2)
IF (x2 .GT. x3) CALL SWAP(x2,x3)
IF (x1 .GT. x2) CALL SWAP(x1,x2)

u = 0.0
IF (x3 .GT. 0.0) u = SQRT(x3)
IF (x2 .LE. 0.0) GO TO 41
IF (x1 .GE. 0.0) GO TO 30
IF (ABS(x1) .GT. x2) GO TO 40
x1 = 0.0

30 x1 = SQRT(x1)
  x2 = SQRT(x2)
  IF (q .GT. 0.0) x1 = -x1
  temp(1) = (( x1 + x2) + u) - b
  temp(2) = ((-x1 - x2) + u) - b
  temp(3) = (( x1 - x2) - u) - b
  temp(4) = ((-x1 + x2) - u) - b
  CALL SELECTSORT(temp)
  IF (ABS(temp(1)) .GE. 0.1*ABS(temp(4))) GO TO 31
  t = temp(2)*temp(3)*temp(4)
  IF (t /= 0.0) temp(1) = e/t

31 z(1) = CMPLX(temp(1), 0.0,DP)
  z(2) = CMPLX(temp(2), 0.0,DP)
  z(3) = CMPLX(temp(3), 0.0,DP)
  z(4) = CMPLX(temp(4), 0.0,DP)
  RETURN

40 v1 = SQRT(ABS(x1))
v2 = 0.0
GO TO 50

41 v1 = SQRT(ABS(x1))
v2 = SQRT(ABS(x2))
IF (q .LT. 0.0) u = -u

50 x = -u - b
y    = v1 - v2
z(1) = CMPLX(x, y,DP)
z(2) = CMPLX(x,-y,DP)
x    = u - b
y    = v1 + v2
z(3) = CMPLX(x, y,DP)
z(4) = CMPLX(x,-y,DP)
RETURN

! The resolvent cubic has complex roots

60 t = DBLE(z(1))
x = 0.0
IF (t .LT. 0.0) THEN
  GO TO 61
ELSE IF (t == 0.0) THEN
  GO TO 70
ELSE
  GO TO 62
END IF
61 h = ABS(DBLE(z(2))) + ABS(AIMAG(z(2)))
IF (ABS(t) .LE. h) GO TO 70
GO TO 80
62 x = SQRT(t)
IF (q .GT. 0.0) x = -x

70 w = SQRT(z(2))
  u  = 2.0*DBLE(w)
  v  = 2.0*ABS(AIMAG(w))
  t  =  x - b
  x1 = t + u
  x2 = t - u
  IF (ABS(x1) .LE. ABS(x2)) GO TO 71
  t  = x1
  x1 = x2
  x2 = t
71 u = -x - b
  h = u*u + v*v
  IF (x1*x1 .LT. 0.01*MIN(x2*x2,h)) x1 = e/(x2*h)
  z(1) = CMPLX(x1, 0.0,DP)
  z(2) = CMPLX(x2, 0.0,DP)
  z(3) = CMPLX(u, v,DP)
  z(4) = CMPLX(u,-v,DP)
  RETURN

80 v = SQRT(ABS(t))
  z(1) = CMPLX(-b, v,DP)
  z(2) = CMPLX(-b,-v,DP)
  z(3) = z(1)
  z(4) = z(2)
  RETURN

!-------------------------------------------------------------------------------!
END SUBROUTINE QuarticRoots
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE SELECTSORT(a)
!-------------------------------------------------------------------------------!
! Reorder the elements of in increasing order.                                  !
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL(KIND=DP),INTENT(IN OUT) :: a(:)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER                      :: j
INTEGER                      :: k(1)
!-------------------------------------------------------------------------------!

!-------------------------------------------------!
! This is a n**2 method.                          !
! It should only be used for small arrays. n<25   !
!-------------------------------------------------!

DO j = 1, SIZE(a)-1
  k = MINLOC(a(j:))
  IF (j .NE. k(1)) CALL SWAP(a(k(1)),a(j))
END DO
RETURN

!-------------------------------------------------------------------------------!
END SUBROUTINE SELECTSORT
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE SWAP(a,b)
!-------------------------------------------------------------------------------!
! Interchange the contents of a and b                                           !
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL(KIND=DP),INTENT(INOUT) :: a, b
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL(KIND=DP)               :: t
!-------------------------------------------------------------------------------!

t = b
b = a
a = t
RETURN

!-------------------------------------------------------------------------------!
END SUBROUTINE SWAP
!===============================================================================!
!
!
!
!-------------------------------------------------------------------------------!
END MODULE MOD_QuarticRoots
!===============================================================================!
