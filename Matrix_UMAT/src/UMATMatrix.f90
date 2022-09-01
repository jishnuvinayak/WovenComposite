!DEC$ FREEFORM
!===============================================================================
! Material routine for matrix of woven composite.
!===============================================================================
include 'ABQinterface.f90'
include 'typeMaterialProperties_Matrix.f90'
include 'typeStateVariables_Matrix.f90'
include 'UMATMatrix_Utilities.f90'
SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT,&
                DRPLDE,DRPLDT,&
                STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,&
                NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,&
                CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)

    USE ABQINTERFACE
    USE Type_State_variables_matrix
    USE Type_material_prop_matrix
    USE Matrix_utilities , only:compute_stress,&
							  perform_viscoleastic_correction,compute_yield_function,&
							  compute_viscoelastic_material_tangent,perform_viscoplastic_correction,&
							  compute_viscoplastic_material_tangent
    IMPLICIT NONE
    !===============================================================================
    !UMAT interface variable declarations
    !===============================================================================
    !--------------------Variables to be defined------------------------------------
    real(kind=AbqRK)::STRESS(NTENS),&       !Array of stress
                      DDSDDE(NTENS,NTENS),& !Material tangent matrix
                      STATEV(NSTATV)        !Array of solution dependent internal state variables
    real(kind=AbqRK)::SSE,SPD,SCD           ! Specific elastic strain energy, plastic and creep dissipation
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
                         
    real(kind=AbqRK)::STRAN(NTENS),&    	!Array with strain at the begining of time inc
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
    ! Material parameter variables
    Type(material_prop_matrix)::mat_properties
    !History variables @ n
    Type(state_variables_matrix)::state_var_previous
    !Trial state variables
    Type(state_variables_matrix)::state_var_trial
    !Updated variables @ n+1
    real(kind=AbqRK)::Strain_total_n2(NTENS)
    Type(state_variables_matrix)::state_var_current
	!Tangent variables
    real(kind=AbqRK)::C_t(NTENS,NTENS)	!Material tangent
    !Misc variables
	real(kind=AbqRK)::yield_f			!yield criteria value
	integer(kind=AbqIK)::i,j,num_NR_iterations
    !Parameters
    real(kind=AbqRK),PARAMETER:: allowed_ZeroTolerance	=	1e-6_AbqRK
    
    !===========================================================================
    !UMAT logic -> Start
    !===========================================================================
    
    !Get the material properties from PROPS
    !===========================================================================
    !E_e -> PROPS[1]
    !Nu - > PROPS[2]
    !Nnum Kelvin branches(N) -> PROPS[3]
    !Viscoelastic_prop -> PROPS[3:2+N*2]
    !ViscoPlasticDamageProps[R0,K,n,H,m,S,Beta] -> PROPS[3+N*2: 9+N*2]					   
    !
    !===========================================================================
    call mat_properties%allocation(PROPS,NTENS)
    
    !Fetch solution dependent state variables
    !===========================================================================
    !Strain_Vp    -> STATEV[1:NTENS]
    !Strain_Ve    -> STATEV[NTENS+1:NTENS*(N+1)]
    !Plastic_mult -> STATEV[NTENS*(N+1)+1]
    !Damage_var   -> STATEV[NTENS*(N+1)+2]
    !===========================================================================
    call state_var_previous%allocation(NTENS,mat_properties%N_k)
   
    state_var_previous%Strain_Vp 									= STATEV(1:NTENS)
    forall(i=1:mat_properties%N_k)state_var_previous%Strain_Ve(i,:)	= STATEV((NTENS*i)+1:NTENS*(i+1))
    state_var_previous%Plastic_mult 								= STATEV(NTENS*(mat_properties%N_k+1)+1)
    state_var_previous%Damage_var 									= STATEV(NTENS*(mat_properties%N_k+1)+2)
    !LineSearchFlag													= STATEV(NSTATV-1)
    
!    write(6,*),'-------------In Time step',TIME(2),'-------------  -------------'
    !Variable initializations
    PNEWDT					= 1.0_AbqRK
    C_t 					= 0.0_AbqRK
    num_NR_iterations 		= 0.0_AbqIK
    
    !Allocate and initialize current and trial state variables
    call state_var_trial%allocation(NTENS,mat_properties%N_k)
    call state_var_current%allocation(NTENS,mat_properties%N_k)
    
    !Strain @ n+1
    Strain_total_n2 = STRAN+DSTRAN
   
    !-----------------------Step 1- Viscoelastic correction---------------------
    ! (assumption : no plastic/damage evolution)
    
    !Perfrom viscoleastic correction : compute new trial values of Strain_Ve_n @ n+1 ie Strain_Ve_trial_n2 
    call perform_viscoleastic_correction(Strain_total_n2,state_var_previous,state_var_trial,mat_properties,mat_properties%N_k,NTENS,DTIME)
											
    
    !Compute updates stress (Stress_trial_n2) using corrected viscoelastic strains(Strain_Ve_n2)
    call compute_stress(mat_properties,Strain_total_n2,state_var_trial%Strain_Ve,state_var_previous%Strain_Vp,state_var_previous%Damage_var,state_var_trial%Stress,NTENS)
   
    !Calculate the yield function value(yield_f)
    yield_f = compute_yield_function(state_var_trial,state_var_previous,mat_properties,NTENS,NDI)
	
    !check yield condition
    if(yield_f <= allowed_ZeroTolerance) then
		!No yielding in this step. Accept the corrected stress state and state variables
		state_var_current%Strain_Ve 	= state_var_trial%Strain_Ve
		state_var_current%Strain_Vp 	= state_var_trial%Strain_Vp
		state_var_current%Plastic_mult	= state_var_trial%Plastic_mult
		state_var_current%Damage_var 	= state_var_trial%Damage_var
		state_var_current%Stress 		= state_var_trial%Stress
		
		!Compute the material tangent(C_t)
		call compute_viscoelastic_material_tangent(mat_properties,state_var_current,C_t,DTIME,mat_properties%N_k,NTENS)
		PNEWDT					= 1.3_AbqRK
		
	
    else
	!------------------Step 3- Viscoelastic/plastic correction-----------------------------
		
		!Perfrom viscoplastic correction : compute new values of state_var_trial
		call perform_viscoplastic_correction(Strain_total_n2,state_var_previous,state_var_trial,mat_properties,yield_f,DTIME,PNEWDT,NDI,mat_properties%N_k,NTENS,num_NR_iterations)
		
		!Compute updates stress (state_var_trial%Stress) using corrected trail states
		call compute_stress(mat_properties,Strain_total_n2,state_var_trial%Strain_Ve,state_var_trial%Strain_Vp,state_var_trial%Damage_var,state_var_trial%Stress,NTENS)
																			
		!update corrected state
		state_var_current%Strain_Ve 	= state_var_trial%Strain_Ve
		state_var_current%Strain_Vp 	= state_var_trial%Strain_Vp
		state_var_current%Plastic_mult	= state_var_trial%Plastic_mult
		state_var_current%Damage_var 	= state_var_trial%Damage_var
		state_var_current%Stress 		= state_var_trial%Stress
		
		!Compute the material tangent(C_t)
		call compute_viscoplastic_material_tangent(Strain_total_n2,state_var_current,state_var_previous,mat_properties,C_t,NDI,DTIME,mat_properties%N_k,NTENS)
										
    end if
    
    !update stress and material tangent
    STRESS = state_var_current%Stress
    DDSDDE = C_t 
    
	!Update the solution dependent state variables
	STATEV(1:NTENS)   							  					= state_var_current%Strain_Vp(:)
	forall(i=1:mat_properties%N_k) STATEV((NTENS*i)+1:NTENS*(i+1)) 	= state_var_current%Strain_Ve(i,:)
	STATEV(NTENS*(mat_properties%N_k+1)+1) 						  	= state_var_current%Plastic_mult
	STATEV(NTENS*(mat_properties%N_k+1)+2) 						  	= state_var_current%Damage_var
	STATEV(NSTATV) 						  							= num_NR_iterations
   
    !Deallocate allocatable memory
    call state_var_previous%deallocation()
    call state_var_current%deallocation()
    call state_var_trial%deallocation()
    call mat_properties%deallocation()
    
    END SUBROUTINE UMAT
