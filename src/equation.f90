!===============================================================================!
MODULE MOD_Equation
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
PRIVATE
!-------------------------------------------------------------------------------!
INTERFACE ExactFunction
  MODULE PROCEDURE ExactFunction
END INTERFACE

INTERFACE SourceTerms
  MODULE PROCEDURE SourceTerms
END INTERFACE

INTERFACE BoundaryConditions
  MODULE PROCEDURE BoundaryConditions
END INTERFACE

INTERFACE TimeStep
  MODULE PROCEDURE TimeStep
END INTERFACE

INTERFACE RiemannSolver
  MODULE PROCEDURE RiemannSolver
END INTERFACE

INTERFACE EvaluateFlux1D
  MODULE PROCEDURE EvaluateFlux1D
END INTERFACE

INTERFACE ConsToPrim
  MODULE PROCEDURE ConsToPrim
END INTERFACE

INTERFACE PrimToCons
  MODULE PROCEDURE PrimToCons
END INTERFACE
!-------------------------------------------------------------------------------!
PUBLIC :: ExactFunction
PUBLIC :: SourceTerms
PUBLIC :: BoundaryConditions
PUBLIC :: TimeStep
PUBLIC :: RiemannSolver
PUBLIC :: EvaluateFlux1D
PUBLIC :: ConsToPrim
PUBLIC :: PrimToCons
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
SUBROUTINE ExactFunction(WhichInitialCondition,t,x,Cons)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: PI
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X0
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X1
USE MOD_FiniteVolume2D_vars,ONLY: MESH_SX
USE MOD_FiniteVolume2D_vars,ONLY: Kappa
USE MOD_FiniteVolume2D_vars,ONLY: KappaM1
USE MOD_FiniteVolume2D_vars,ONLY: KappaP1
USE MOD_FiniteVolume2D_vars,ONLY: sKappaM1
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState1
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState2
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState3
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState4
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
INTEGER,INTENT(IN) :: WhichInitialCondition
REAL,INTENT(IN)    :: t
REAL,INTENT(IN)    :: x(1:nDims)
REAL,INTENT(OUT)   :: Cons(1:nVar)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL                :: Prim(1:nVar)
REAL                :: Prim_in(1:nVar), Prim_out(1:nVar), Prim_am(1:nVar)
REAL                :: delta_rho, delta_vx, delta_vy, xm, ym
REAL                :: r, xc(2), x0, xt, r0, y0
REAL                :: rho0, p0, vx0, vy0, nx, ny, om, phi
REAL                :: sigma, a_shear, v_shear, am, rho1
CHARACTER(LEN=255)  :: ErrorMessage
!-------------------------------------------------------------------------------!

Cons = 0.0
Prim = 0.0
SELECT CASE(WhichInitialCondition)
  !----------------------------------------------------------------------!
  ! [200] Constant State                                                 !
  !----------------------------------------------------------------------!
  CASE(200)
    Prim(1:nVar) = PrimRefState1(1:nVar)

    CALL PrimToCons(Prim,Cons)
  !----------------------------------------------------------------------!
  ! [211] Riemann Problem                                                !
  !----------------------------------------------------------------------!
  CASE(211)
    xm = MESH_X0(1)+0.5*MESH_SX(1)
    ym = MESH_X0(2)+0.5*MESH_SX(2)
    IF      ((x(1) .GE. xm) .AND. (x(2) .GE. ym)) THEN
      Prim(1:nVar) = PrimRefState1(1:nVar)
    ELSE IF ((x(1) .LT. xm) .AND. (x(2) .GE. ym)) THEN
      Prim(1:nVar) = PrimRefState2(1:nVar)
    ELSE IF ((x(1) .LT. xm) .AND. (x(2) .LT. ym)) THEN
      Prim(1:nVar) = PrimRefState3(1:nVar)
    ELSE IF ((x(1) .GE. xm) .AND. (x(2) .LT. ym)) THEN
      Prim(1:nVar) = PrimRefState4(1:nVar)
    END IF

    CALL PrimToCons(Prim,Cons)
  !----------------------------------------------------------------------!
  ! [214] Double Mach Reflection Problem                                 !
  !----------------------------------------------------------------------!
  CASE(214)
    Prim_in(1:nVar)  = PrimRefState1(1:nVar)
    Prim_out(1:nVar) = PrimRefState2(1:nVar)
    x0 = 1.0/6.0
    xt = x0 + (1.0/SQRT(3.0))*x(2)
    IF (x(1) .LT. xt) THEN
      Prim(1:nVar) = Prim_in(1:nVar)
    ELSE
      Prim(1:nVar) = Prim_out(1:nVar)
    END IF

    CALL PrimToCons(Prim,Cons)
  !----------------------------------------------------------------------!
  ! [216] Cloud-Shock Interaction                                        !
  !----------------------------------------------------------------------!
  CASE(216)
    IF (x(1) .LE. MESH_X0(1)+0.8*MESH_SX(1)) THEN
      Prim(1:nVar) = PrimRefState1(1:nVar)
    ELSE
      Prim(1:nVar) = PrimRefState2(1:nVar)
    END IF
    xm    = MESH_X0(1)+0.7*MESH_SX(1)
    ym    = MESH_X0(2)+0.5*MESH_SX(2)
    xc(1) = x(1)-xm
    xc(2) = x(2)-ym
    r     = SQRT(xc(1)**2 + xc(2)**2)
    IF (r .LE. 0.15) THEN
      Prim(1:nVar) = PrimRefState3(1:nVar)
    END IF
    CALL PrimToCons(Prim,Cons)
  !----------------------------------------------------------------------!
  ! [217] Kelvin-Helmholtz Instability                                   !
  !----------------------------------------------------------------------!
  CASE(217)
    y0        = 0.50 !0.25
    sigma     = 0.10
    a_shear   = 0.01
    v_shear   = 0.50
    am        = 0.10
    rho0      = 0.505
    rho1      = 0.495

    IF (x(2) .GT. 0.0) THEN
      delta_rho = rho0 + rho1*TANH((x(2)-y0)/a_shear)
      delta_vx  = +v_shear*TANH((x(2)-y0)/a_shear)
      delta_vy  = +am*v_shear*SIN(2.0*PI*x(1))*EXP(-(x(2)-y0)**2/sigma**2)
    ELSEIF (x(2) .LE. 0.0) THEN
      delta_rho = rho0 - rho1*TANH((x(2)+y0)/a_shear)
      delta_vx  = -v_shear*TANH((x(2)+y0)/a_shear)
      delta_vy  = -am*v_shear*SIN(2.0*PI*x(1))*EXP(-(x(2)+y0)**2/sigma**2)
    END IF

    Prim(1) = delta_rho
    Prim(2) = delta_vx
    Prim(3) = delta_vy
    Prim(4) = 1.0

    CALL PrimToCons(Prim,Cons)
  CASE DEFAULT
    ErrorMessage = "Exact function not specified"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!-------------------------------------------------------------------------------!
END SUBROUTINE ExactFunction
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE SourceTerms(t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: S
USE MOD_FiniteVolume2D_vars,ONLY: MeshBary
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER          :: ii, jj
!-------------------------------------------------------------------------------!

S = 0.0

!-------------------------------------------------------------------------------!
END SUBROUTINE SourceTerms
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE BoundaryConditions(t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X0
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X1
USE MOD_FiniteVolume2D_vars,ONLY: MeshBary
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState1
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState2
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState3
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState4
USE MOD_FiniteVolume2D_vars,ONLY: BoundaryConditionsType
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: V
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)    :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj
INTEGER            :: idx_vx, idx_vy
REAL               :: x0, xc, xt
REAL               :: Prim_in(1:nVar), Prim_out(1:nVar)
REAL               :: Cons_in(1:nVar), Cons_out(1:nVar)
REAL               :: ConsRefState1(1:nVar), ConsRefState2(1:nVar)
CHARACTER(LEN=255) :: ErrorMessage
!-------------------------------------------------------------------------------!

idx_vx = 2
idx_vy = 3

!------------------------------!
! Left Boundary Conditions     !
!------------------------------!
SELECT CASE(BoundaryConditionsType(4))
  CASE(1) ! Periodic
    DO jj=1,nElemsY
      DO ii=1,nGhosts
        U(1:nVar,-nGhosts+ii,jj) = U(1:nVar,nElemsX-nGhosts+ii,jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO jj=1,nElemsY
      DO ii=1,nGhosts
        U(1:nVar,-nGhosts+ii,jj) = U(1:nVar,nGhosts-ii+1,jj)
      END DO
    END DO
  CASE(3) ! Inflow
    Prim_in(1:nVar) = PrimRefState4(1:nVar)
    CALL PrimToCons(Prim_in(1:nVar),Cons_in(1:nVar))
    DO jj=1,nElemsY
      DO ii=1,nGhosts
        U(1:nVar,-nGhosts+ii,jj) = Cons_in(1:nVar)
      END DO
    END DO
  CASE(4) ! Outflow
    Prim_out(1:nVar) = PrimRefState4(1:nVar)
    CALL PrimToCons(Prim_out(1:nVar),Cons_out(1:nVar))
    DO jj=1,nElemsY
      DO ii=1,nGhosts
        U(1:nVar,-nGhosts+ii,jj) = Cons_out(1:nVar)
      END DO
    END DO
  CASE(5) ! Reflecting
    DO jj=1,nElemsY
      DO ii=1,nGhosts
        U(1:nVar,-nGhosts+ii,jj) = U(1:nVar,nGhosts-ii+1,jj)
        U(idx_vx,-nGhosts+ii,jj) =-U(idx_vx,nGhosts-ii+1,jj)
      END DO
    END DO
  CASE(11) ! Double Mach Reflection
    ErrorMessage = "Boundary condition defined only for faces 1 and 3"
    WRITE(*,*) ErrorMessage
    STOP
  CASE DEFAULT
    ErrorMessage = "Boundary condition not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!------------------------------!
! Right Boundary Conditions    !
!------------------------------!
SELECT CASE(BoundaryConditionsType(2))
  CASE(1) ! Periodic
    DO jj=1,nElemsY
      DO ii=1,nGhosts
        U(1:nVar,nElemsX+ii,jj) = U(1:nVar,ii,jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO jj=1,nElemsY    
      DO ii=1,nGhosts
        U(1:nVar,nElemsX+ii,jj) = U(1:nVar,nElemsX-ii+1,jj)
      END DO
    END DO
  CASE(3) ! Inflow
    Prim_in(1:nVar) = PrimRefState2(1:nVar)
    CALL PrimToCons(Prim_in(1:nVar),Cons_in(1:nVar))
    DO jj=1,nElemsY
      DO ii=1,nGhosts
        U(1:nVar,nElemsX+ii,jj) = Cons_in(1:nVar)
      END DO
    END DO
  CASE(4) ! Outflow
    Prim_out(1:nVar) = PrimRefState2(1:nVar)
    CALL PrimToCons(Prim_out(1:nVar),Cons_out(1:nVar))
    DO jj=1,nElemsY
      DO ii=1,nGhosts
        U(1:nVar,nElemsX+ii,jj) = Cons_out(1:nVar)
      END DO
    END DO
  CASE(5) ! Reflecting
    DO jj=1,nElemsY
      DO ii=1,nGhosts
        U(1:nVar,nElemsX+ii,jj) = U(1:nVar,nElemsX-ii+1,jj)
        U(idx_vx,nElemsX+ii,jj) =-U(idx_vx,nElemsX-ii+1,jj)
      END DO
    END DO
  CASE(11) ! Double Mach Reflection
    ErrorMessage = "Boundary condition defined only for faces 1 and 3"
    WRITE(*,*) ErrorMessage
    STOP
  CASE DEFAULT
    ErrorMessage = "Boundary condition not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!------------------------------!
! Top Boundary Conditions      !
!------------------------------!
SELECT CASE(BoundaryConditionsType(3))
  CASE(1) ! Periodic
    DO ii=1,nElemsX
      DO jj=1,nGhosts
        U(1:nVar,ii,nElemsY+jj) = U(1:nVar,ii,jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO ii=1,nElemsX
      DO jj=1,nGhosts
        U(1:nVar,ii,nElemsY+jj) = U(1:nVar,ii,nElemsY-jj+1)
      END DO
    END DO
  CASE(3) ! Inflow
    Prim_in(1:nVar) = PrimRefState3(1:nVar)
    CALL PrimToCons(Prim_in(1:nVar),Cons_in(1:nVar))
    DO ii=1,nElemsX
      DO jj=1,nGhosts
        U(1:nVar,ii,nElemsY+jj) = Cons_in(1:nVar)
      END DO
    END DO
  CASE(4) ! Outflow
    Prim_out(1:nVar) = PrimRefState3(1:nVar)
    CALL PrimToCons(Prim_out(1:nVar),Cons_out(1:nVar))
    DO ii=1,nElemsX
      DO jj=1,nGhosts
        U(1:nVar,ii,nElemsY+jj) = Cons_out(1:nVar)
      END DO
    END DO
  CASE(5) ! Reflecting
    DO ii=1,nElemsX
      DO jj=1,nGhosts
        U(1:nVar,ii,nElemsY+jj) = U(1:nVar,ii,nElemsY-jj+1)
        U(idx_vy,ii,nElemsY+jj) =-U(idx_vy,ii,nElemsY-jj+1)
      END DO
    END DO
  CASE(11) ! Double Mach Reflection
    Prim_in(1:nVar)  = PrimRefState1(1:nVar)
    Prim_out(1:nVar) = PrimRefState3(1:nVar)
    CALL PrimToCons(Prim_in(1:nVar),Cons_in(1:nVar))
    CALL PrimToCons(Prim_out(1:nVar),Cons_out(1:nVar))
    DO ii=1,nElemsX
      DO jj=1,nGhosts
        x0 = 1.0/6.0
        xc = MeshBary(1,ii,nElemsY)
        xt = x0 + (1.0+2.0*0.4984*t)/SQRT(3.0)
        IF (xc .LT. xt) THEN
          U(1:nVar,ii,nElemsY+jj) = Cons_in(1:nVar)
        ELSE
          U(1:nVar,ii,nElemsY+jj) = Cons_out(1:nVar)
        END IF
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!------------------------------!
! Bottom Boundary Conditions   !
!------------------------------!
SELECT CASE(BoundaryConditionsType(1))
  CASE(1) ! Periodic
    DO ii=1,nElemsX
      DO jj=1,nGhosts
        U(1:nVar,ii,-nGhosts+jj) = U(1:nVar,ii,nElemsY-nGhosts+jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO ii=1,nElemsX
      DO jj=1,nGhosts
        U(1:nVar,ii,-nGhosts+jj) = U(1:nVar,ii,nGhosts-jj+1)
      END DO
    END DO
  CASE(3) ! Inflow
    Prim_in(1:nVar) = PrimRefState1(1:nVar)
    CALL PrimToCons(Prim_in(1:nVar),Cons_in(1:nVar))
    DO ii=1,nElemsX
      DO jj=1,nGhosts
        U(1:nVar,ii,-nGhosts+jj) = Cons_in(1:nVar)
      END DO
    END DO
  CASE(4) ! Outflow
    Prim_out(1:nVar) = PrimRefState1(1:nVar)
    CALL PrimToCons(Prim_out(1:nVar),Cons_out(1:nVar))
    DO ii=1,nElemsX
      DO jj=1,nGhosts
        U(1:nVar,ii,-nGhosts+jj) = Cons_out(1:nVar)
      END DO
    END DO
  CASE(5) ! Reflecting
    DO ii=1,nElemsX
      DO jj=1,nGhosts
        U(1:nVar,ii,-nGhosts+jj) = U(1:nVar,ii,nGhosts-jj+1)
        U(idx_vy,ii,-nGhosts+jj) =-U(idx_vy,ii,nGhosts-jj+1)
      END DO
    END DO
  CASE(11) ! Double Mach Reflection
    Prim_in(1:nVar) = PrimRefState1(1:nVar)
    CALL PrimToCons(Prim_in(1:nVar),Cons_in(1:nVar))
    DO ii=1,nElemsX
      DO jj=1,nGhosts
        x0 = 1.0/6.0
        xc = MeshBary(1,ii,1)
        IF (xc .LT. x0) THEN
          U(1:nVar,ii,-nGhosts+jj) = Cons_in(1:nVar)
        ELSE
          U(1:nVar,ii,-nGhosts+jj) = U(1:nVar,ii,nGhosts-jj+1)
          U(idx_vy,ii,-nGhosts+jj) =-U(idx_vy,ii,nGhosts-jj+1)
        END IF
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

DO jj=-nGhosts,nElemsY+nGhosts+1
  DO ii=-nGhosts,nElemsX+nGhosts+1
    CALL ConsToPrim(U(1:nVar,ii,jj),V(1:nVar,ii,jj))
  END DO
END DO

!-------------------------------------------------------------------------------!
END SUBROUTINE BoundaryConditions
!===============================================================================!
!
!
!
!===============================================================================!
FUNCTION TimeStep()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: CFL
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: LambdaMaxX
USE MOD_FiniteVolume2D_vars,ONLY: LambdaMaxY
USE MOD_FiniteVolume2D_vars,ONLY: MIN_TIMESTEP
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL    :: TimeStep
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL    :: FastestWaveX, FastestWaveY
REAL    :: Prim(1:nVar)
INTEGER :: ii, jj
!-------------------------------------------------------------------------------!

LambdaMaxX = 0.0
LambdaMaxY = 0.0
TimeStep = HUGE(1.0)

DO jj=1,nElemsY
  DO ii=1,nElemsX
    CALL ConsToPrim(U(1:nVar,ii,jj),Prim(1:nVar))
    CALL WaveSpeeds2D(Prim(1:nVar),FastestWaveX,FastestWaveY)
    LambdaMaxX = MAX(LambdaMaxX,ABS(FastestWaveX))
    LambdaMaxY = MAX(LambdaMaxY,ABS(FastestWaveY))
    TimeStep  = MIN(TimeStep,MESH_DX(1)/LambdaMaxX,MESH_DX(2)/LambdaMaxY)
  END DO
END DO

TimeStep = CFL*TimeStep

IF (TimeStep .LT. MIN_TIMESTEP) THEN
  TimeStep = MIN_TIMESTEP
END IF

!-------------------------------------------------------------------------------!
END FUNCTION TimeStep
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE WaveSpeeds1D(Prim,slowest,fastest)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)           :: Prim(1:nVar)
REAL,INTENT(OUT),OPTIONAL :: slowest
REAL,INTENT(OUT),OPTIONAL :: fastest
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL                      :: rho, vx, vy, p
REAL                      :: c, v2, c2, qx
!-------------------------------------------------------------------------------!

rho = Prim(1)
vx  = Prim(2)
vy  = Prim(3)
p   = Prim(4)

c  = EOS_SoundSpeed(rho,p)
c2 = c*c
v2 = vx*vx + vy*vy
qx = (1.0-v2)*(1.0-vx*vx-c2*(v2-vx*vx))

IF (PRESENT(slowest)) THEN
  slowest = vx*(1.0-c2)-c*SQRT(qx)
  slowest = slowest/(1.0-v2*c2)
END IF
IF (PRESENT(fastest)) THEN
  fastest = vx*(1.0-c2)+c*SQRT(qx)
  fastest = fastest/(1.0-v2*c2)
END IF

!-------------------------------------------------------------------------------!
END SUBROUTINE WaveSpeeds1D
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE WaveSpeeds2D(Prim,fastestx,fastesty)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Prim(1:nVar)
REAL,INTENT(OUT) :: fastestx
REAL,INTENT(OUT) :: fastesty
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: rho, vx, vy, p
REAL             :: c, v2, c2, qx, qy
!-------------------------------------------------------------------------------!

rho = Prim(1)
vx  = Prim(2)
vy  = Prim(3)
p   = Prim(4)

c  = EOS_SoundSpeed(rho,p)
c2 = c*c
v2 = vx*vx + vy*vy
qx = (1.0-v2)*(1.0-vx*vx-c2*(v2-vx*vx))
qy = (1.0-v2)*(1.0-vy*vy-c2*(v2-vy*vy))

fastestx = vx*(1.0-c2)+c*SQRT(qx)
fastestx = fastestx/(1.0-v2*c2)

fastesty = vy*(1.0-c2)+c*SQRT(qy)
fastesty = fastesty/(1.0-v2*c2)

!-------------------------------------------------------------------------------!
END SUBROUTINE WaveSpeeds2D
!===============================================================================!
!
!
!
!===============================================================================!
FUNCTION EOS_SoundSpeed(rho,p)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: Kappa
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: rho, p
REAL            :: EOS_SoundSpeed
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL               :: h
CHARACTER(LEN=255) :: ErrorMessage
!-------------------------------------------------------------------------------!

h = EOS_Enthalpy(rho,p)

EOS_SoundSpeed = (Kappa*p)/(rho*h)

IF (EOS_SoundSpeed .LT. 0.0) THEN
  ErrorMessage = "Negative speed of sound"
  WRITE(*,*) ErrorMessage
  STOP
END IF

EOS_SoundSpeed = SQRT(EOS_SoundSpeed)

!-------------------------------------------------------------------------------!
END FUNCTION EOS_SoundSpeed
!===============================================================================!
!
!
!
!===============================================================================!
FUNCTION EOS_Enthalpy(rho,p)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: Kappa,KappaM1
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: rho, p
REAL            :: EOS_Enthalpy
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!

EOS_Enthalpy = 1.0 + (Kappa/KappaM1)*(p/rho)

!-------------------------------------------------------------------------------!
END FUNCTION EOS_Enthalpy
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE ConsToPrim(Cons, Prim)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: MAXITER
USE MOD_FiniteVolume2D_vars,ONLY: Kappa, sKappaM1
USE MOD_FiniteVolume2D_vars,ONLY: EPS, ACCURACY
USE MOD_FiniteVolume2D_vars,ONLY: MIN_DENSITY, MIN_ENERGY, MIN_MOMENTUM
USE MOD_FiniteVolume2D_vars,ONLY: MIN_PRESSURE, MIN_SPEED, MAX_SPEED
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Cons(1:nVar)
REAL,INTENT(OUT) :: Prim(1:nVar)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL               :: D, Sx, Sy, E
REAL               :: rho, vx, vy, p
REAL               :: S, v, v2
REAL               :: LF, LF2
REAL               :: h, temp
REAL               :: fp, dfp
LOGICAL            :: STATUS
INTEGER            :: iter
CHARACTER(LEN=255) :: ErrorMessage
!-------------------------------------------------------------------------------!

D  = Cons(1)
Sx = Cons(2)
Sy = Cons(3)
E  = Cons(4)

S  = SQRT(Sx*Sx + Sy*Sy)

!------------------------------!
! Unphysical values detector   !
!------------------------------!
IF (D .LT. MIN_DENSITY) THEN
  D = MIN_DENSITY
END IF
IF (ABS(Sx) .LT. MIN_MOMENTUM) THEN
  Sx = 0.0
END IF
IF (ABS(Sy) .LT. MIN_MOMENTUM) THEN
  Sy = 0.0
END IF
IF (E .LT. MIN_ENERGY) THEN
  E = SQRT(D*D + S*S + EPS)
END IF
temp = SQRT(D*D + S*S)
IF (E .LT. temp) THEN
  E = SQRT(temp*temp + EPS)
END IF

p = S - E
p = MAX(ABS(p),MIN_PRESSURE)

iter   = 0
STATUS = .FALSE.
DO WHILE (.NOT. STATUS)
  v   = S/(E + p)
  LF  = 1.0-v*v
  LF  = MAX(LF, EPS)
  LF  = 1.0/LF
  LF  = SQRT(LF)
  rho = D/LF
  h   = EOS_Enthalpy(rho, p)

  fp  = D*h*LF - E - p
  dfp = (Kappa*sKappaM1)*LF*LF &
      - (S*S*LF*LF*LF/((E+p)*(E+p)*(E+p)))*(D+(2.0*Kappa*sKappaM1)*p*LF)-1.0
  p   = p - fp/dfp

  IF(ABS(fp) .LT. ACCURACY) THEN
    STATUS = .TRUE.
  END IF

  iter = iter + 1
  IF (iter .EQ. MAXITER) THEN
    CALL ConsToPrim_ExactFormula(Cons,Prim)
    RETURN
  END IF
END DO

IF (ABS(v) .LT. MIN_SPEED) THEN
  vx = 0.0
  vy = 0.0
ELSE
  vx = v*Sx/S
  vy = v*Sy/S
END IF
IF (ABS(v) .GT. MAX_SPEED) THEN
  vx = vx*MAX_SPEED/v
  vy = vy*MAX_SPEED/v
  v  = vx*vx + vy*vy
  v  = SQRT(v)
END IF

IF (p .NE. p) THEN
  ErrorMessage = "NaN found while recovering pressure in ConsToPrim"
  WRITE(*,*) ErrorMessage
  STOP
END IF

IF (p .LT. 0.0) THEN
  ErrorMessage = "Negative pressure in ConsToPrim"
  WRITE(*,*) ErrorMessage
  STOP
END IF

Prim(1) = D/LF
Prim(2) = vx
Prim(3) = vy
Prim(4) = p

!-------------------------------------------------------------------------------!
END SUBROUTINE ConsToPrim
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE ConsToPrim_ExactFormula(Cons,Prim)
!-------------------------------------------------------------------------------!
! Based on ryu2006a:                                                            !
!   Dongsu Ryu, Indranil Chattopadhyay and Eunwoo Choi                          !
!   Equation of state in numerical relativistic hydrodynamics                   !
!   The Astrophysical Journal Supplement Series, 166:410Y420, 2006              !
!-------------------------------------------------------------------------------!
USE MOD_QuarticRoots,       ONLY: QuarticRoots
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: MAXITER
USE MOD_FiniteVolume2D_vars,ONLY: Kappa, KappaM1
USE MOD_FiniteVolume2D_vars,ONLY: EPS, ACCURACY
USE MOD_FiniteVolume2D_vars,ONLY: MIN_DENSITY, MIN_ENERGY, MIN_MOMENTUM
USE MOD_FiniteVolume2D_vars,ONLY: MIN_PRESSURE, MIN_SPEED, MAX_SPEED
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Cons(1:nVar)
REAL,INTENT(OUT) :: Prim(1:nVar)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: D, Sx, Sy, E
REAL             :: rho, vx, vy, p
REAL             :: S, v, v2
REAL             :: LF, LF2
REAL             :: temp
REAL             :: aa(0:4)
COMPLEX          :: ss(1:4)
!-------------------------------------------------------------------------------!

D  = Cons(1)
Sx = Cons(2)
Sy = Cons(3)
E  = Cons(4)

S  = SQRT(Sx*Sx + Sy*Sy)

!------------------------------!
! Unphysical values detector   !
!------------------------------!
IF (D .LT. MIN_DENSITY) THEN
  D = MIN_DENSITY
END IF
IF (E .LT. MIN_ENERGY) THEN
  E = SQRT(D*D + S*S + EPS)
END IF
temp = SQRT(D*D + S*S)
IF (E .LT. temp) THEN
  E = SQRT(temp*temp + EPS)
END IF

!------------------------------!
! Polynomial coefficients      !
!------------------------------!
temp = (D*D + S*S)*KappaM1*KappaM1
temp = 1.0/temp

aa(0) = S*S*temp
aa(1) = -2.0*Kappa*S*E*temp
aa(2) = (Kappa*Kappa*E*E + 2.0*KappaM1*S*S - KappaM1*KappaM1*D*D)*temp
aa(3) = -2.0*Kappa*KappaM1*S*E*temp
aa(4) = 1.0

!------------------------------!
! The root finder              !
!------------------------------!
CALL QuarticRoots(aa,ss)
v  = ss(2)
v2 = v*v
IF (ABS(v) .LT. MIN_SPEED) THEN
  vx = 0.0
  vy = 0.0
ELSE
  vx = v*Sx/S
  vy = v*Sy/S
END IF
IF (v .GT. MAX_SPEED) THEN
  v  = 1.0/v
  vx = vx*v*MAX_SPEED
  vy = vy*v*MAX_SPEED
  v2 = vx*vx + vy*vy
  v  = SQRT(v2)
END IF

LF  = 1.0/SQRT(1.0-v2)
rho = D/LF
p   = KappaM1*((E - Sx*vx -Sy*vy)-rho)

IF (p .LT. MIN_PRESSURE) THEN
  p = MIN_PRESSURE
END IF

!------------------------------------------------------------!
! Variables Transformation                                   !
!------------------------------------------------------------!
Prim(1) = rho
Prim(2) = vx
Prim(3) = vy
Prim(4) = p

!-------------------------------------------------------------------------------!
END SUBROUTINE ConsToPrim_ExactFormula
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE PrimToCons(Prim, Cons)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: MIN_RHO, MIN_PRESSURE, MIN_SPEED, MAX_SPEED
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Prim(1:nVar)
REAL,INTENT(OUT) :: Cons(1:nVar)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: rho, vx, vy, p
REAL             :: v, v2, LF, LF2, h
!-------------------------------------------------------------------------------!

rho = Prim(1)
vx  = Prim(2)
vy  = Prim(3)
p   = Prim(4)

IF (ABS(vx) .LT. MIN_SPEED) THEN
  vx = 0.0
END IF
IF (ABS(vy) .LT. MIN_SPEED) THEN
  vy = 0.0
END IF

v2 = vx*vx + vy*vy
v  = SQRT(v2)
!------------------------------!
! Unphysical Values Corrector  !
!------------------------------!
IF (rho .LT. MIN_RHO) THEN
  rho = MIN_RHO
END IF

IF (v .GT. MAX_SPEED) THEN
  v  = 1.0/v
  vx = vx*v*MAX_SPEED
  vy = vy*v*MAX_SPEED
  v2 = vx*vx + vy*vy
  v  = SQRT(v2)
END IF

IF (p .LT. MIN_PRESSURE) THEN
  p = MIN_PRESSURE
END IF

!------------------------------!
! The Variables Mapping        !
!------------------------------!
h   = EOS_Enthalpy(rho,p)
LF2 = 1.0/(1.0-v2)
LF  = SQRT(LF2)

Cons(1) = rho*LF
Cons(2) = rho*h*LF2*vx
Cons(3) = rho*h*LF2*vy
Cons(4) = rho*h*LF2-p

!-------------------------------------------------------------------------------!
END SUBROUTINE PrimToCons
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE EvaluateFlux1D(Prim,Flux)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: Kappa, sKappaM1
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Prim(1:nVar)
REAL,INTENT(OUT) :: Flux(1:nVar)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: rho, vx, vy, p
REAL             :: v2, LF2, LF, h
!-------------------------------------------------------------------------------!

rho = Prim(1)
vx  = Prim(2)
vy  = Prim(3)
p   = Prim(4)

v2  = vx*vx + vy*vy
h   = EOS_Enthalpy(rho,p)
LF2 = 1.0/(1.0-v2)
LF  = SQRT(LF2)

Flux(1) = rho*LF*vx
Flux(2) = rho*h*LF2*vx*vx + p
Flux(3) = rho*h*LF2*vx*vy
Flux(4) = rho*h*LF2*vx

!-------------------------------------------------------------------------------!
END SUBROUTINE EvaluateFlux1D
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE RiemannSolver(PrimL,PrimR,NormVect,TangVect,Flux)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: nDims
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: PrimL(1:nVar,1:nGPs)
REAL,INTENT(IN)  :: PrimR(1:nVar,1:nGPs)
REAL,INTENT(IN)  :: NormVect(1:nDims,1:nGPs)
REAL,INTENT(IN)  :: TangVect(1:nDims,1:nGPs)
REAL,INTENT(OUT) :: Flux(1:nVar,1:nGPs)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: PrimLL(1:nVar,1:nGPs), PrimRR(1:nVar,1:nGPs)
REAL             :: ConsLL(1:nVar,1:nGPs), ConsRR(1:nVar,1:nGPs)
INTEGER          :: iGP
!-------------------------------------------------------------------------------!

DO iGP=1,nGPs
  ! Rotating the vector quantities       !
  PrimLL(1,iGP) = PrimL(1,iGP)
  PrimLL(2,iGP) = NormVect(1,iGP)*PrimL(2,iGP) + NormVect(2,iGP)*PrimL(3,iGP)
  PrimLL(3,iGP) = TangVect(1,iGP)*PrimL(2,iGP) + TangVect(2,iGP)*PrimL(3,iGP)
  PrimLL(4,iGP) = PrimL(4,iGP)

  PrimRR(1,iGP) = PrimR(1,iGP)
  PrimRR(2,iGP) = NormVect(1,iGP)*PrimR(2,iGP) + NormVect(2,iGP)*PrimR(3,iGP)
  PrimRR(3,iGP) = TangVect(1,iGP)*PrimR(2,iGP) + TangVect(2,iGP)*PrimR(3,iGP)
  PrimRR(4,iGP) = PrimR(4,iGP)

  CALL PrimToCons(PrimLL(1:nVar,iGP),ConsLL(1:nVar,iGP))
  CALL PrimToCons(PrimRR(1:nVar,iGP),ConsRR(1:nVar,iGP))  

  CALL RiemannSolverByRusanov(&
    ConsLL(1:nVar,iGP),ConsRR(1:nVar,iGP),&
    PrimLL(1:nVar,iGP),PrimRR(1:nVar,iGP),Flux(1:nVar,iGP))

  ! Rotating back the vector quantities  !
  Flux(2:3,iGP) = NormVect(1:nDims,iGP)*Flux(2,iGP) &
                + TangVect(1:nDims,iGP)*Flux(3,iGP)
END DO

!-------------------------------------------------------------------------------!
END SUBROUTINE RiemannSolver
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE RiemannSolverByRusanov(ConsL,ConsR,PrimL,PrimR,Flux)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: ConsL(1:nVar), ConsR(1:nVar)
REAL,INTENT(IN)  :: PrimL(1:nVar), PrimR(1:nVar)
REAL,INTENT(OUT) :: Flux(1:nVar)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: FluxL(1:nVar), FluxR(1:nVar)
REAL             :: LambdaMax, fastestL, fastestR
!-------------------------------------------------------------------------------!

CALL EvaluateFlux1D(PrimL,FluxL)
CALL EvaluateFlux1D(PrimR,FluxR)
CALL WaveSpeeds1D(PrimL,fastest=fastestL)
CALL WaveSpeeds1D(PrimR,fastest=fastestR)

LambdaMax = MAX(ABS(fastestL),ABS(fastestR))

Flux = 0.5*((FluxL + FluxR) - LambdaMax*(ConsR - ConsL))

!-------------------------------------------------------------------------------!
END SUBROUTINE RiemannSolverByRusanov
!===============================================================================!
!
!
!
!===============================================================================!
END MODULE MOD_Equation
!-------------------------------------------------------------------------------!
