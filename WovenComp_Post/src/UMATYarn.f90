!DEC$ FREEFORM
!===============================================================================
! Material routine for yarn of woven composite.
!===============================================================================
!include 'ABQinterface.f90'
!include 'UMATYarn_Utilities.f90'
SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT,&
                DRPLDE,DRPLDT,&
                STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,&
                NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,&
                CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)

    USE ABQINTERFACE
    USE Yarn_utilities , only:assemble_initial_stiffness,assemble_eshelby_tensor,&
							  compute_interaction_tensor,compute_stiffness_reduction,&
							  compute_damage_activation_criteria,compute_residual_phi,&
							  correct_trial_state, compute_material_tangent
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
    real(kind=AbqRK)::C_0(NTENS,NTENS)						!Initial material stifness in Viogt notation
    real(kind=AbqRK)::S_E(NTENS,NTENS)						!Eshelby tensor in Viogt notation
    real(kind=AbqRK)::R(2) 									!Initial threshold in transverse tension and in plane shear
    real(kind=AbqRK)::WB_P(2) 								!Weibull law parameters
    real(kind=AbqRK)::gamma_c_infinity 						!Micro crack saturation parameter
    real(kind=AbqRK)::Inelasticity_param(2) 				!Inelasticity parameters ( a22,a12)
    ! History variables @ n
    real(kind=AbqRK)::Strain_inelastic_n(NTENS) 			!Inelastic strain from previous time step n
    real(kind=AbqRK)::gamma_c_n 							!Micro crack density from previous time step n
    real(kind=AbqRK)::Hc_sup								!Supremum value of stress criterion calculated in the entire load history
    !Trial state variables
    real(kind=AbqRK)::TrialStress(NTENS) 					!Trial stress for stress prediction
    real(kind=AbqRK)::TrialStrain_inelastic(NTENS)			!Trial inelastic strain.
    real(kind=AbqRK)::Trialgamma_c 							!Trial Micro crack density
    !Corrected variables
    real(kind=AbqRK)::CorrectedStress(NTENS) 				!Corrected stress for stress prediction
    real(kind=AbqRK)::CorrectedStrain_inelastic(NTENS)		!Corrected inelastic strain.
    real(kind=AbqRK)::Correctedgamma_c						!Corrected Micro crack density
    !Updated variables @ n+1
    real(kind=AbqRK)::Strain_inelastic_n2(NTENS)			!Updated inelastic strain.
    real(kind=AbqRK)::gamma_c_n2							!Updated Micro crack density
    real(kind=AbqRK)::Strain_total_n2(NTENS) 
	! Tangent variables
    real(kind=AbqRK)::D_stiffness_reduction(NTENS,NTENS)	!Stiffness reduction tensor
    real(kind=AbqRK)::C_t(NTENS,NTENS)						!Material tangent
    !Misc variables
    real(kind=AbqRK)::H_c									!Damage activation criterion for current load step
    real(kind=AbqRK)::phi_gamma 							!Residual for gama
    real(kind=AbqRK)::T_c(NTENS,NTENS) 						!interaction tensor
    !Parameters
    real(kind=AbqRK),PARAMETER:: allowed_ZeroTolerance	=	0.0000001_AbqRK
    
    !===========================================================================
    !UMAT logic -> Start
    !===========================================================================
    !Variable initializations
    D_stiffness_reduction	= 0.0_AbqRK
    C_0 					= 0.0_AbqRK
    S_E 					= 0.0_AbqRK
    T_c 					= 0.0_AbqRK
    H_c 					= 0.0_AbqRK
    C_t 					= 0.0_AbqRK
    PNEWDT					= 1.5_AbqRK
    
    !Strain @ n+1
    Strain_total_n2 = STRAN+DSTRAN
    
    !Get the material properties from PROPS
    !===========================================================================
    !C_0 -> PROPS[1:5]  order : C1111,C1122,C2222,C2233,C1212 (transverse isotropy)
    !S_E -> PROPS[6:17] order :SE1111,SE1122,SE1133,SE2211,SE2222,SE2233,SE3311,
    !						   SE3322,SE3333,SE1212,SE1313,SE2323( 12 non zero components)
    !R	  -> PROPS[18:19] order : R22,R12
    !WB_P -> PROPS[20:21] order : S,Beta
    !gamma_c_infinity -> PROPS[22]
    !Inelasticity_param -> PROPS[23:24] order : a22,a12
    !===========================================================================
    call assemble_initial_stiffness(C_0,PROPS(1:5),NTENS)
    call assemble_eshelby_tensor(S_E,PROPS(6:17),NTENS)
    R(1:2)=PROPS(18:19)
    WB_P(1:2)=PROPS(20:21)
    gamma_c_infinity = PROPS(22)
    Inelasticity_param(1:2) = PROPS(23:24)
    
    !Fetch solution dependent state variables
    !===========================================================================
    !Strain_inelastic_n -> STATEV[1:NTENS]
    !gama_c_n -> STATEV[NTENS+1]
    !Hc_sup -> STATEV[NSTATV]
    !===========================================================================
    Strain_inelastic_n(1:NTENS) = STATEV(1:NTENS)
    gamma_c_n = STATEV(NTENS+1)
    Hc_sup = STATEV(NSTATV)
    
    !Compute interaction tensor once -> T_c
    call compute_interaction_tensor(S_E,T_c,NTENS)
    !Compute stiffness reductionn tesnor D(gama_c) -> D_stiffness_reduction
    call compute_stiffness_reduction(gamma_c_n,D_stiffness_reduction,C_0,T_c,NTENS)
   
    
    !-----------------------Stress Prediction-----------------------------------
    
    !Calculate the trial state
    TrialStress = reshape(MATMUL((C_0 - D_stiffness_reduction),&
						reshape((Strain_total_n2-Strain_inelastic_n),&
									(/NTENS,1/))),(/NTENS/))
    TrialStrain_inelastic = Strain_inelastic_n
    Trialgamma_c = gamma_c_n
    !Compute damage activation criterion using the trial stress state -> H_c
    
    call compute_damage_activation_criteria(H_c,TrialStress,Trialgamma_c,R)
    
    !Check activation of damage
    
    if(H_c < 1.0_AbqRK) then
		!Damage is not yet activated. In elastic region. Accept trial state
		STRESS = TrialStress
		DDSDDE = C_0-D_stiffness_reduction
		Strain_inelastic_n2 = TrialStrain_inelastic
		gamma_c_n2 = Trialgamma_c
		Hc_sup = H_c
		
    else
		! Damage is active and need to check if it evolves or not.
		
		!Compute the residual for damage with the trial values -> phi_gamma
		call compute_residual_phi(H_c,gamma_c_n,phi_gamma,WB_P,gamma_c_infinity)
		
		if( phi_gamma <= allowed_ZeroTolerance) then
			! No evolution of gama. Can accept the trial state
			STRESS = TrialStress
			DDSDDE = C_0-D_stiffness_reduction
			Strain_inelastic_n2 = TrialStrain_inelastic
			gamma_c_n2 = Trialgamma_c
		else
	!------------------------Stress correction----------------------------------
			!Need to compute the updated value of gamma, inealstic strain and correct the trial state					 
			call correct_trial_state(TrialStrain_inelastic,Trialgamma_c,TrialStress,&
									 CorrectedStrain_inelastic,Correctedgamma_c,&
									 CorrectedStress, H_c,WB_p,R,gamma_c_infinity,&
									 Inelasticity_param,Strain_total_n2,C_0,T_c,NTENS,PNEWDT)
			
			!update the stress and internal state variables with corrected state
			STRESS = CorrectedStress
			Strain_inelastic_n2 = CorrectedStrain_inelastic
			gamma_c_n2 = Correctedgamma_c
			
			!Compute the material tangent -> C_t
			call compute_material_tangent(gamma_c_n2,STRESS,C_0,T_c,&
										Strain_total_n2,Strain_inelastic_n2,&
										Inelasticity_param,R,WB_P,C_t,&
										gamma_c_infinity,NTENS)
			
			!update the material stiffness
			DDSDDE = C_t
			
		end if
    end if
    
    !Update the solution dependent state variables
    
    STATEV(1:NTENS) = Strain_inelastic_n2(1:NTENS)
    STATEV(NTENS+1) = gamma_c_n2
    STATEV(NTENS+2:NTENS+8) = STRESS(:)
    if (H_c > Hc_sup) STATEV(NSTATV) = H_c
    
    
    END SUBROUTINE UMAT
