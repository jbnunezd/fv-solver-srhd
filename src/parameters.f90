!===============================================================================!
MODULE MOD_Parameters
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
PRIVATE
!-------------------------------------------------------------------------------!
INTERFACE InitializeParameters
  MODULE PROCEDURE InitializeParameters
END INTERFACE
!-------------------------------------------------------------------------------!
PUBLIC :: InitializeParameters
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
SUBROUTINE InitializeParameters()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: PI
USE MOD_FiniteVolume2D_vars,ONLY: CFL
USE MOD_FiniteVolume2D_vars,ONLY: TEnd
USE MOD_FiniteVolume2D_vars,ONLY: Kappa
USE MOD_FiniteVolume2D_vars,ONLY: KappaM1
USE MOD_FiniteVolume2D_vars,ONLY: KappaP1
USE MOD_FiniteVolume2D_vars,ONLY: sKappaM1
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X0
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X1
USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction
USE MOD_FiniteVolume2D_vars,ONLY: ReconstructionFix
USE MOD_FiniteVolume2D_vars,ONLY: WhichOutput
USE MOD_FiniteVolume2D_vars,ONLY: nOutputFiles
USE MOD_FiniteVolume2D_vars,ONLY: InitialCondition
USE MOD_FiniteVolume2D_vars,ONLY: BoundaryConditionsType
USE MOD_FiniteVolume2D_vars,ONLY: VarNameVisu
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState1
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState2
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState3
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState4
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
CHARACTER(LEN=255) :: ErrorMessage
!-------------------------------------------------------------------------------!

InitialCondition = 211

PrimRefState1 = (/1.0,0.0,0.0,3.6/)
PrimRefState2 = (/1.0,0.0,0.0,3.6/)
PrimRefState3 = (/1.0,0.0,0.0,3.6/)
PrimRefState4 = (/1.0,0.0,0.0,3.6/)

SELECT CASE(InitialCondition)
  CASE(200) ! Constant State
    TEnd    = 1.0
    Kappa   = 1.4
    nElemsX = 100
    nElemsY = 100
    MESH_X0 = (/0.0,0.0/)
    MESH_X1 = (/1.0,1.0/)
    BoundaryConditionsType = (/1,1,1,1/)
  CASE(211) ! Riemann Problem
    TEnd    = 1.0
    Kappa = 1.66666666666667
    nElemsX = 200
    nElemsY = 200
    MESH_X0 = (/-1.0,-1.0/)
    MESH_X1 = (/+1.0,+1.0/)
    PrimRefState1 = (/0.035145216124503,0.00,0.00,0.162931056509027/)
    PrimRefState2 = (/0.10,0.70,0.00,1.00/)
    PrimRefState3 = (/0.50,0.00,0.00,1.00/)
    PrimRefState4 = (/0.10,0.00,0.70,1.00/)
    BoundaryConditionsType = (/2,2,2,2/)
  CASE(214) ! Double Mach Reflection Problem
    TEnd    = 4.0
    Kappa   = 1.4
    nElemsX = 400
    nElemsY = 100
    MESH_X0 = (/0.0,0.0/)
    MESH_X1 = (/4.0,1.0/)
    PrimRefState1 = (/8.564,0.4247*SIN(PI/3.0),-0.4247*COS(PI/3.0),0.3808/)
    PrimRefState2 = (/1.4,0.0,0.0,0.0025/)
    PrimRefState3 = (/1.4,0.0,0.0,0.0025/)
    PrimRefState4 = (/8.564,0.4247*SIN(PI/3.0),-0.4247*COS(PI/3.0),0.3808/)
    BoundaryConditionsType = (/11,4,11,3/)
  CASE(216) ! Cloud-Shock Interaction
    TEnd    = 4.00
    Kappa   = 1.66666666666667
    nElemsX = 800
    nElemsY = 400
    MESH_X0 = (/0.0,0.0/)
    MESH_X1 = (/2.0,1.0/)
    PrimRefState1 = (/1.00,0.00,0.00,0.05/)
    PrimRefState2 = (/1.865225080631180,-0.196781107378299,0.00,0.15/)
    PrimRefState3 = (/3.1538,0.00,0.00,0.05/)
    BoundaryConditionsType = (/2,2,2,2/)
  CASE(217) ! Kelvin-Helmholtz Instability
    TEnd    = 5.00
    Kappa   = 1.33333333333333
    nElemsX = 400
    nElemsY = 800
    MESH_X0 = (/-0.5,-1.0/)
    MESH_X1 = (/+0.5,+1.0/)
    BoundaryConditionsType = (/1,1,1,1/)
  CASE DEFAULT
    ErrorMessage = "Initial condition not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

CFL      = 0.95
KappaM1  = Kappa-1.0
KappaP1  = Kappa+1.0
sKappaM1 = 1.0/KappaM1

Reconstruction    = 3
ReconstructionFix = 2

WhichOutput  = 1
nOutputFiles = 100

VarNameVisu(1) = 'Density'
VarNameVisu(2) = 'VelocityX'
VarNameVisu(3) = 'VelocityY'
VarNameVisu(4) = 'Pressure'

!-------------------------------------------------------------------------------!
END SUBROUTINE InitializeParameters
!===============================================================================!
!
!
!
!===============================================================================!
END MODULE MOD_Parameters
!-------------------------------------------------------------------------------!
