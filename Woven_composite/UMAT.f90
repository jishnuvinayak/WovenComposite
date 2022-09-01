!DEC$ FREEFORM
!===============================================================================
! UMAT like interface for switching between material routines of yarn and matrix
! of woven composites.
!===============================================================================
SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT,&
                DRPLDE,DRPLDT,&
                STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,&
                NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,&
                CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)

    USE ABQINTERFACE
    USE Yarn, UMATYarn=>UMAT        		!Module with yarn material routine
    USE Matrix, UMATMatrix=>UMAT    		!Module with matrix material routine
    IMPLICIT NONE
    !===============================================================================
    !UMAT interface variable declarations
    !===============================================================================
    !--------------------Variables to be defined------------------------------------
    real(kind=AbqRK)::STRESS(NTENS),&       !Array of stress
                      DDSDDE(NTENS,NTENS),& !Material tangent matrix
                      STATEV(NSTATV)        !Array of solution dependent internal state variables
    real(kind=AbqRK)::SSE,SPD,SCD           !Specific elastic strain energy, plastic and creep dissipation
    real(kind=AbqRK)::RPL                   !Volumetric heat generation per unit time
    real(kind=AbqRK)::DDSDDT(NTENS),&       !Variation of stress increment wrt tempertature
                      DRPLDE(NTENS),&       !Variation of RPL wrt to strain increment
                      DRPLDT                !Variation of RPL wrt to temperature

    !--------------------Variables that can be updated------------------------------
    real(kind=AbqRK)::PNEWDT                !Ratio of suggested time increment to the increment currently used

    !--------------------Variables passed for information---------------------------
    integer(kind=AbqIK)::NOEL,NPT,&         !Element number, Integration point number
                         KSPT,KINC,LAYER    !Section number, Increment number, Layer number
    integer(kind=AbqIK)::JSTEP(4)           !Array with step details
    integer(kind=AbqIK)::NDI,NSHR,NTENS,&   !Num of direct stress component,shear stress component, total stress components
                         NSTATV,&           !Number of solution dependent internal state variables
                         NPROPS             !Number of materail properties
                         
    real(kind=AbqRK),intent(in)::STRAN(NTENS),&    	!Array with strain at the begining of time inc
                      DSTRAN(NTENS)     	!Array with strain increments
    real(kind=AbqRK)::TIME(2),&         	!Value of step time and total time at begning of current inc
                      DTIME             	!Time increment
    real(kind=AbqRK)::TEMP,DTEMP        	!Temperature and temperature increment at start of current time inc
    real(kind=AbqRK)::PREDEF(1),&       	!Array of interpolated values of predefined field variables.
                      DPRED(1)          	!Array of increments of predefined field variables.		  
    real(kind=AbqRK)::PROPS(NPROPS)     	!Array with user defined material properties
    real(kind=AbqRK)::COORDS(3),&       	!Array with co-ordinates of the current point
                      DROT(3,3),&       	!Rotation increment matrix
                      CELENT            	!Characteristic element length
    real(kind=AbqRK)::DFGRD0(3,3),&     	!Array with deformation gradient at the begining of inc
                      DFGRD1(3,3)       	!Array with deformation gradient at the end of inc
    character *80 CMNAME                	!User defined material name
    !===============================================================================
    !Local variable declaration.
    !===============================================================================
    integer(kind=AbqRK)::element_material,i 	!To store material reference passed in from micro UEL
    real(kind=AbqRK)::e_yarn_1(3),e_yarn_2(3) 	!Yarn orientations ( e1 and e2)
    real(kind=AbqRK)::rotation_matrix(3,3),identity(3,3)	!Rotation matrix
    real(kind=AbqRK)::STRESS_R(NTENS) 		!Rotated stress
    real(kind=AbqRK)::STRAN_R(NTENS),&
					  DSTRAN_R(NTENS),DDSDDE_R(NTENS,NTENS) 		!Rotated strain
    
    !===============================================================================
    !Subroutine logic -> Start
    !===============================================================================
    !PNEWDT					= 25.0_AbqRK
        
    ! Check the material ref passed by micro UEL to select the material routine
    element_material = int(PROPS(1))
    if (element_material== 1) then
    
        !Call matrix material routine if the material point corresponds to matrix
        call UMATMatrix(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT,&
                         DRPLDE,DRPLDT,STRAN,DSTRAN,TIME,&
                         DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,&
                         NDI,NSHR,NTENS,NSTATV,[PROPS(2:NPROPS)],NPROPS-1,COORDS,&
                         DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,&
                         LAYER,KSPT,JSTEP(1),KINC)
    else if (element_material == 2) then
		!Calcuate the rotation matrix from yarn orientation
		e_yarn_1 = PROPS(2:4)
		e_yarn_2 = PROPS(5:7)
		rotation_matrix = compute_rotation_matrix(e_yarn_1,e_yarn_2)
		
		!Rotate stress
		call ROTSIG(STRESS,rotation_matrix,STRESS_R,1,NDI,NSHR)
		
		!Rotate strain and dstrain
		call ROTSIG(STRAN,rotation_matrix,STRAN_R,2,NDI,NSHR)
		call ROTSIG(DSTRAN,rotation_matrix,DSTRAN_R,2,NDI,NSHR) 
		    
        !Call yarn material routine if the material point corresponds to yarn
        call UMATYarn(STRESS_R,STATEV,DDSDDE_R,SSE,SPD,SCD,RPL,DDSDDT,&
                        DRPLDE,DRPLDT,STRAN_R,DSTRAN_R,TIME,&
                        DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,&
                        NDI,NSHR,NTENS,NSTATV,PROPS(8:NPROPS),NPROPS-7,COORDS,&
                        DROT,PNEWDT,CELENT,&
                        DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
        
        !Rotate the stress back
		call ROTSIG(STRESS_R,transpose(rotation_matrix),STRESS,1,NDI,NSHR)
		!Rotate material tangent to GCS
        DDSDDE =rotate_fourth_order_tensor(DDSDDE_R,rotation_matrix,NTENS)
		
    else
        !Throw error message for undefined materials
        call STDB_ABQERR(-3,"Error: Material not defined", 0, 0.0, " ")
    end if

END SUBROUTINE UMAT

