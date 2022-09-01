!DEC$ FREEFORM
!===============================================================================
! Seperate utility module for yarn UMAT. Conatins all helper subroutines and
! functions called in UMAT yarn.
!===============================================================================
MODULE Yarn_utilities
	USE ABQINTERFACE
    IMPLICIT NONE
    Private
	Public :: assemble_initial_stiffness, 			assemble_eshelby_tensor,&
			  I_4,compute_interaction_tensor,		compute_stiffness_reduction,&
			  compute_damage_activation_criteria,	compute_residual_phi,&
			  compute_A_gamma_Hc,					compute_lambda_inelastic_strain,&
			  compute_derivative_A0_gamma, 			compute_B_Hc_gamma,correct_trial_state,&
			  compute_derivative_stress_strain, 	compute_derivative_Hc_strain,&
			  compute_derivative_D_gamma, 			compute_derivative_stress_gamma,&
			  compute_material_tangent, 			invert_matrix
	Contains
	!-----------------------Procedures--------------------------------------
	
	Subroutine assemble_initial_stiffness(C_0,PROPS,NTENS)
		!=======================================================================
		!Subroutine to assembel intial material stiffness as 6*6 matrix from 
		!PROPS
		!	Inputs:-
		!	PROPS	-> Non zero components of initial material stiffness
		!            order : C1111,C1122,C2222,C2233,C1212
		!    
		!	Outputs:-
		!    C_0		-> Assembled 6x6 material stiffness in voigt form
		!=======================================================================
		IMPLICIT NONE
		!Dummy variable declarations
		real(kind=AbqRK),intent(in)::PROPS(:)
		integer(kind=AbqIK),intent(in)::NTENS
		real(kind=AbqRK),intent(out)::C_0(NTENS,NTENS)
		
		!Variable initialization
		C_0 = 0.0_AbqRK
		
		!For tranversely isotropic material ( all non zero components)
		!1 row
		C_0(1,1) = PROPS(1)
		C_0(1,2) = PROPS(2)
		C_0(1,3) = C_0(1,2)
		!2nd row
		C_0(2,1) = C_0(1,2)
		C_0(2,2) = PROPS(3)
		C_0(2,3) = PROPS(4)
		! 3rd row
		C_0(3,1) = C_0(1,3)
		C_0(3,2) = C_0(2,3)
		C_0(3,3) = C_0(2,2)
		! 4th row
		C_0(4,4) = PROPS(5)
		! 5th row
		C_0(5,5) = C_0(4,4)
		! 6th row
		C_0(6,6) = (0.5_AbqRK)*(C_0(2,2)-C_0(2,3))
		
		
	End Subroutine
	
	Subroutine assemble_eshelby_tensor(S_E,PROPS,NTENS)
		!=======================================================================
		!Subroutine to assembel eshelby tensor as 6*6 matrix from PROPS
		!
		!	Inputs:-
		!	PROPS	-> Non zero components of initial material stiffness
		!
		!            order :SE1111,SE1122,SE1133,SE2211,SE2222,SE2233,SE3311,
		!			 SE3322,SE3333,SE1212,SE1313,SE2323( 12 non zero components)
		!    
		!	Outputs:-
		!    S_E	-> Assembled 6x6 material stiffness in voigt form
		!=======================================================================
		IMPLICIT NONE
		!Dummy variable declarations
		real(kind=AbqRK),intent(in)::PROPS(:)
		integer(kind=AbqIK),intent(in)::NTENS
		real(kind = AbqRK),intent(out)::S_E(NTENS,NTENS)
		
		!Variable initialization
		S_E = 0.0_AbqRK
		
		!For tranversely isotropic material ( all non zero components)
		!1 row
		S_E(1,1) = PROPS(1)
		S_E(1,2) = PROPS(2)
		S_E(1,3) = PROPS(3)
		!2ns row
		S_E(2,1) = PROPS(4)
		S_E(2,2) = PROPS(5)
		S_E(2,3) = PROPS(6)
		! 3rd row
		S_E(3,1) = PROPS(7)
		S_E(3,2) = PROPS(8)
		S_E(3,3) = PROPS(9)
		! 4th row
		S_E(4,4) = PROPS(10)
		! 5th row
		S_E(5,5) = PROPS(11)
		! 6th row
		S_E(6,6) = PROPS(12)
			
	End Subroutine
	
	Subroutine compute_interaction_tensor(S_E,T_c,NTENS)
		!=======================================================================
		!Subroutine to interaction tensor from S_E
		!
		!	Inputs:-
		!	S_E		-> Eshelby tensor
		!	NTENS 	-> Number of stress compoenents
		!    
		!	Outputs:-
		!   T_c		-> Interaction tensor
		!=======================================================================
		IMPLICIT NONE
		!Dummy variable declarations
		integer(kind=AbqIK),intent(in)::NTENS
		real(kind=AbqRK),intent(in)::S_E(:,:)
		real(kind=AbqRK),intent(out)::T_c(:,:)
		
		!Local variables
		real(kind=AbqRK)::TempHolderT_c(NTENS,NTENS) 	!Temp holder for interaction tensor
		real(kind=AbqRK)::Identity_tensor(NTENS,NTENS)	!Fourth order identity tensor in voigt form
	
		!Variable initializations
		T_c = 0.0_AbqRK
		TempHolderT_c = 0.0_AbqRK
		Identity_tensor = 0.0_AbqRK
		
		!-------------Compute interaction tensor--------------------------------
		Identity_tensor = I_4()
		TempHolderT_c = Identity_tensor - S_E
		!T_c is the inverse of TempHolderT_c
		call invert_matrix(TempHolderT_c,T_c)
			
	End Subroutine
	
	Subroutine compute_stiffness_reduction(gamma_c_n,D_stiffness_reduction,C_0,T_c,NTENS)
		!=======================================================================
		!Subroutine to compute stiffness reduction tensor D(gama)
		!
		!	Inputs:-
		!	gamma_c_n	-> micro crack density @ time step n
		!	C_0			-> Initial material stiffness
		!	T_c			-> Interaction tensor
		!	NTENS 		-> Number of stress compoenents
		!
		!	Outputs:-
		!   D_stiffness_reduction	-> Stiffness reduction tensor
		!=======================================================================
		IMPLICIT NONE
		!Dummy variable declarations
		integer(kind=AbqIK),intent(in)::NTENS
		real(kind=AbqRK),intent(in)::gamma_c_n
		real(kind=AbqRK),intent(in)::C_0(:,:)
		real(kind=AbqRK),intent(in)::T_c(:,:)
		real(kind=AbqRK),intent(out)::D_stiffness_reduction(:,:)
		
		!Local variables
		real(kind=AbqRK)::A_0(NTENS,NTENS) 			 	!Strain localisation tensor
		real(kind=AbqRK)::TempHolderD_0(NTENS,NTENS) 	!Temp holder for matrix multiplication
		!Variable initializations
		D_stiffness_reduction 	= 0.0_AbqRK
		A_0 					= 0.0_AbqRK
		TempHolderD_0 			= 0.0_AbqRK
		
		!--------Compute strain localisation tensor of virgin medium------------
		A_0	= compute_A0_gamma(gamma_c_n,T_c,NTENS)
		
		!---------Compute the stiffness reduction tensor------------------------
		! D_0 = gamma*C_0:T_c:A_0
		TempHolderD_0			= matmul(C_0,T_c)
		D_stiffness_reduction	= gamma_c_n*matmul(TempHolderD_0,A_0)
		
	End Subroutine
	
	Subroutine compute_damage_activation_criteria(H_c,TrialStress,Trialgamma_c,R)	
		!=======================================================================
		!Subroutine to compute the damage activation criteria H_c
		!
		!	Inputs:-
		!	Trialgamma_c	-> Trial micro crack density
		!	TrialStress 	-> Trial stress
		!	R 				-> Material threshold in 22 and 11 direction
		!
		!	Outputs:-
		!   H_c	-> Damage activation criterion H_c
		!=======================================================================
		IMPLICIT NONE
		!Dummy variable declarations
		real(kind=AbqRK),intent(in)::Trialgamma_c
		real(kind=AbqRK),intent(in)::TrialStress(:)
		real(kind=AbqRK),intent(in)::R(2)
		real(kind=AbqRK),intent(out)::H_c
		
		!Local variables
		real(kind=AbqRK)::LocalStress(2)
		
		!variable initializations
		H_c= 0.0_AbqRK
		LocalStress = 0.0_AbqRK
		
		!Calculate the local stress components S22 and S12
		LocalStress(1) = (1.0_AbqRK/(1.0_AbqRK-Trialgamma_c))*TrialStress(2)
		LocalStress(2) = (1.0_AbqRK/(1.0_AbqRK-Trialgamma_c))*TrialStress(4)
		
		!Compute H_c
		H_c = sqrt((LocalStress(1)/R(1))**2+(LocalStress(2)/R(2))**2)
		
		
	End Subroutine
	
	Subroutine compute_residual_phi(H_c,gamma_c_n,phi_gamma,WB_P,gamma_c_infinity)
		!=======================================================================
		!Subroutine to compute the residual phi_gama, to check for 
		!evolution of damage.
		!
		!	Inputs:-
		!	gama_c_n			-> micro crack density
		!	H_c 				-> yield criteria
		!	WB_P 				-> Weibul law parameters
		!	gamma_c_infinity 	-> micro-crack saturation parameter
		!
		!	Outputs:-
		!   phi_gama			-> residual phi_gama
		!=======================================================================
		IMPLICIT NONE
		!Dummy variable declarations
		real(kind=AbqRK),intent(in)::gamma_c_n
		real(kind=AbqRK),intent(in)::H_c,gamma_c_infinity
		real(kind=AbqRK),intent(in)::WB_P(2)
		real(kind=AbqRK),intent(out)::phi_gamma
		
		!local variable declaration
		real(kind=AbqRK)::exponent_value,H_c_temp
		
		!Variable initializations
		phi_gamma 		= 0.0_AbqRK
		exponent_value	= 0.0_AbqRK
		H_c_temp  = H_c 
		!Condition for macually brackets
		if (H_c <= 1.0_AbqRK) H_c_temp = 1.0_AbqRK
		
		!Calculate Phi_gamma
		exponent_value = exp(-((H_c_temp-1)/WB_P(1))**WB_P(2))
		phi_gamma = (gamma_c_infinity* (1-exponent_value))-gamma_c_n
		
	End Subroutine
	
	Subroutine correct_trial_state(TrialStrain_inelastic,Trialgamma_c,TrialStress,&
									 CorrectedStrain_inelastic,Correctedgamma_c,&
									 CorrectedStress, H_c,WB_p,R,gamma_c_infinity,&
									 Inelasticity_param,Strain_total_n2,C_0,T_c,NTENS,PNEWDT)
		!=======================================================================
		!Subroutine to correct the trial state using cutting plane algorithm
		!	Inputs:-
		!	TrialStrain_inelastic		-> Trial inelastic strain
		!	Trialgamma_c				-> Trial micro crack density
		!	TrialStress					-> Trial stress
		!	H_c							-> Damage activation criteria
		!	WB_p 						-> Weibull law parameters
		!	R							-> Threshold material parameters
		!	gamma_c_infinity			-> Saturation micro crack density parameter
		!	Strain_total_n2 			-> Total strain @ time step n+1
		!	C_0							-> Initial material stiffness
		!	T_c							-> Interaction tensor
		!	NTENS						-> Number of stress compoenents
		!	PNEWDT						-> Ratio of suggested time increment to the increment currently used
		!	
		!
		!	Outputs:-
		!	CorrectedStrain_inelastic	-> Corrected inelastic strain
		!	Correctedgamma_c			-> Corrected micro crack density
		!	CorrectedStress				-> Corrected stress
		!=======================================================================
		IMPLICIT NONE
		!Dummy variable declarations
		integer(kind=AbqIK),intent(in)::NTENS
		real(kind=AbqRK),intent(in)::Trialgamma_c
		real(kind=AbqRK),intent(in)::TrialStrain_inelastic(:)
		real(kind=AbqRK),intent(in)::TrialStress(:)
		real(kind=AbqRK),intent(in)::Strain_total_n2(:)
		real(kind=AbqRK),intent(in)::C_0(:,:)
		real(kind=AbqRK),intent(in)::H_c,gamma_c_infinity
		real(kind=AbqRK),intent(in)::WB_P(:)
		real(kind=AbqRK),intent(in)::R(:)
		real(kind=AbqRK),intent(in)::Inelasticity_param(:)
		real(kind=AbqRK),intent(in)::T_c(:,:)
		real(kind=AbqRK),intent(inout)::PNEWDT
		real(kind=AbqRK),intent(out)::Correctedgamma_c
		real(kind=AbqRK),intent(out)::CorrectedStrain_inelastic(:)
		real(kind=AbqRK),intent(out)::CorrectedStress(:)
		
		
		
		!Local variable declarations
		real(kind=AbqIK)::index_i 								! Running index for cutting plane loop
		real(kind=AbqRK)::gamma_c_i 							! Micro crack density @ i th itertaion
		real(kind=AbqRK)::phi_gamma_i 							! Value of residual @ i th itertaion
		real(kind=AbqRK)::H_c_i 								! Value of damage activation @ i th itertaion
		real(kind=AbqRK)::delta_gamma_c 						! Delta value of micro crack density
		real(kind=AbqRK)::Strain_inelastic_i(NTENS) 			! Inelastic strain @ i th itertaion
		real(kind=AbqRK)::Stress_i(NTENS) 						! Stress @ i th itertaion
		real(kind=AbqRK)::delta_Strain_inelastic(NTENS) 		! Delta value for inelastic strain
		real(kind=AbqRK)::D_stiffness_reduction(NTENS,NTENS) 	! Stiffness reduction tensor
		real(kind=AbqRK)::A_gamma_Hc 							! Partial derivative of Phi_gamma with Hc
		real(kind=AbqRK)::B_Hc_gamma 							! Partial derivative of Hc with gamma
		real(kind=AbqRK)::A_gamma_gamma 						! Partial derivative of Phi_gamma with gamma
		real(kind=AbqRK)::K_gamma_gamma 						! A scalar value (ref literature for more details)
		real(kind=AbqRK)::Lambda_inelastic(NTENS) 				! From inelastic flow rule
		real(kind=AbqRK),parameter::TOEL = 1e-6_AbqRK 			! Tolerence for convergence check
		integer(kind=AbqIK),parameter::Max_iteration = 25_AbqIK ! Max number of iterations to break the loop
		
		!Variable initializations
		phi_gamma_i 				= 0.0_AbqRK
		Correctedgamma_c			= 0.0_AbqRK
		CorrectedStress				= 0.0_AbqRK
		CorrectedStrain_inelastic	= 0.0_AbqRK
		D_stiffness_reduction 		= 0.0_AbqRK
		A_gamma_Hc					= 0.0_AbqRK
		B_Hc_gamma					= 0.0_AbqRK
		A_gamma_gamma				= 0.0_AbqRK
		K_gamma_gamma				= 0.0_AbqRK
		
		!-----------------Start of cutting plane loop---------------------------
		
		!set intial guess (i =1) = trial state
		gamma_c_i 			= Trialgamma_c
		Strain_inelastic_i 	= TrialStrain_inelastic
		Stress_i 			= TrialStress
		H_c_i 				= H_c
		call compute_residual_phi(H_c_i,Trialgamma_c,phi_gamma_i,WB_P,gamma_c_infinity) !For phi_gamma_i
		do index_i = 1, Max_iteration
			!Compute delta_gamma_c
			delta_gamma_c = 0.0_AbqRK
			
			A_gamma_Hc 		= compute_A_gamma_Hc(gamma_c_infinity,H_c_i,WB_P)
			A_gamma_gamma 	= -1.0_AbqRK
			B_Hc_gamma 		= compute_B_Hc_gamma(H_c_i,C_0,T_c,Stress_i,Strain_total_n2,&
								Strain_inelastic_i,gamma_c_i,R,Inelasticity_param,&
								NTENS)							
			K_gamma_gamma 	= A_gamma_Hc*B_Hc_gamma +A_gamma_gamma
			delta_gamma_c 	= -phi_gamma_i / K_gamma_gamma
			
			!Compute delta_Strain_inelastic
			delta_Strain_inelastic 	= 0.0_AbqRK
			Lambda_inelastic 		= compute_lambda_inelastic_strain(Stress_i,Inelasticity_param,NTENS)
			delta_Strain_inelastic 	= Lambda_inelastic * delta_gamma_c
			
			!Update internal state variables gamma and inelastic strain with delta values
			gamma_c_i 			= gamma_c_i + delta_gamma_c
			Strain_inelastic_i 	= Strain_inelastic_i + delta_Strain_inelastic
			
			!Calculate updated stiffness reduction tensor = D_stiffness_reduction
			call compute_stiffness_reduction(gamma_c_i,D_stiffness_reduction,C_0,T_c,NTENS)
			
			!Calulcated updated stress using updated stiffness reduction and gamma value
			Stress_i = reshape(MATMUL((C_0 - D_stiffness_reduction),&
							reshape(Strain_total_n2-Strain_inelastic_i,&
									(/NTENS,1/))),(/NTENS/))
			
			!Use updated values to re calculate residual = phi_gamma_i
			H_c_i = 0.0_AbqRK
			call compute_damage_activation_criteria(H_c_i,Stress_i,gamma_c_i,R)
			call compute_residual_phi(H_c_i,gamma_c_i,phi_gamma_i,WB_P,gamma_c_infinity)
			
			!Check for convergence
			if (ABS(phi_gamma_i) <= TOEL) then
				Correctedgamma_c 			= gamma_c_i
				CorrectedStress 			= Stress_i
				CorrectedStrain_inelastic 	= Strain_inelastic_i
				!write (*,*), 'Stress correction converged at :',index_i
				if (index_i<=Max_iteration*0.8_AbqRK) PNEWDT=1.5_AbqRK
				EXIT
			end if
			!if not converging report a convergence error and request a new time step
			if (index_i == Max_iteration ) then
				!write (6,*), 'Stress correction not converged'
				PNEWDT=0.5_AbqRK
				Correctedgamma_c 			= Trialgamma_c
				CorrectedStress 			= TrialStress
				CorrectedStrain_inelastic 	= TrialStrain_inelastic
			end if
		end do
		
	End Subroutine
	
	Subroutine compute_material_tangent(gamma_c_updated,stress_updated,C_0,T_c,&
										Strain_total_n2,Strain_inelastic_updated,&
										Inelasticity_param,R,WB_P,C_t,&
										gamma_c_infinity,NTENS)
		!=======================================================================
		!Subroutine to compute the material tangent
		!
		!	Inputs:-
		!	Strain_inelastic_updated	-> Trial inelastic strain
		!	gamma_c_updated				-> Trial micro crack density
		!	stress_updated				-> Trial stress
		!	R							-> Threshold material parameters
		!	WB_P 				-> Weibul law parameters
		!	gamma_c_infinity			-> Saturation micro crack density parameter
		!	Strain_total_n2 			-> Total strain @ time step n+1
		!	C_0							-> Initial material stiffness
		!	T_c							-> Interaction tensor
		!	NTENS						-> Number of stress components
		!	
		!
		!	Outputs:-
		!   C_t							-> Material tangent
		!=======================================================================
		IMPLICIT NONE
		!Dummy variable declarations
		integer(kind=AbqIK),intent(in)::NTENS
		real(kind=AbqRK),intent(in)::gamma_c_updated
		real(kind=AbqRK),intent(in)::stress_updated(:)
		real(kind=AbqRK),intent(in)::C_0(:,:)
		real(kind=AbqRK),intent(in)::T_c(:,:)
		real(kind=AbqRK),intent(in)::R(:)
		real(kind=AbqRK),intent(in)::WB_P(:)
		real(kind=AbqRK),intent(in)::Strain_total_n2(:)
		real(kind=AbqRK),intent(in)::Strain_inelastic_updated(:)
		real(kind=AbqRK),intent(in)::Inelasticity_param(:)
		real(kind=AbqRK),intent(in)::gamma_c_infinity
		real(kind=AbqRK),intent(out)::C_t(:,:)
		
		!Local variable declerations
		real(kind=AbqRK)::B_sigma_epsilon(NTENS,NTENS)  ! Partial derivate of stress with strain
		real(kind=AbqRK)::B_Hc_epsilon(NTENS) 			! Partial derivative of Hc with strain
		real(kind=AbqRK)::B_Hc_gamma					! Partial derivate of Hc with gamma
		real(kind=AbqRK)::A_gamma_Hc 					! Partial derivative of Phi_gamma with Hc
		real(kind=AbqRK)::K_gamma_gamma 				! A scalar value (ref literature for more details)
		real(kind=AbqRK)::A_gamma_gamma 				! Partial derivative of Phi_gamma with gamma
		real(kind=AbqRK)::X_gamma_epsilon(NTENS)		! Ref literature
		real(kind=AbqRK)::B_sigma_gamma(NTENS) 			! Partial derivative of sigma with gamma
		real(kind=AbqRK)::H_c
		!Variable initializations
		B_sigma_epsilon 	= 0.0_AbqRK
		B_Hc_epsilon 		= 0.0_AbqRK
		A_gamma_Hc 			= 0.0_AbqRK
		H_c 				= 0.0_AbqRK
		B_Hc_gamma 			= 0.0_AbqRK
		K_gamma_gamma 		= 0.0_AbqRK
		X_gamma_epsilon 	= 0.0_AbqRK
		B_sigma_gamma 		= 0.0_AbqRK
		
		!Compute the partial derivative of stress by strain ( B_sigma_epsilon)
		B_sigma_epsilon	= compute_derivative_stress_strain(gamma_c_updated,C_0,T_c,NTENS)
		
		!Compute the partial deriavtive of Hc by strain(B_Hc_epsilon)
		B_Hc_epsilon 	= compute_derivative_Hc_strain(gamma_c_updated,&
													stress_updated,C_0,T_c,R,NTENS)
		
		!Compute partial derivative of residual phi with Hc (A_gamma_Hc)
		A_gamma_Hc		= compute_A_gamma_Hc(gamma_c_infinity,H_c,WB_P)
		
		!Compute updated H_c
		call compute_damage_activation_criteria(H_c,stress_updated,gamma_c_updated,R)
		
		!Compute partial derivative of Hc with gamma (B_Hc_gamma)
		B_Hc_gamma 		= compute_B_Hc_gamma(H_c,C_0,T_c,stress_updated,Strain_total_n2,&
								Strain_inelastic_updated,gamma_c_updated,R,Inelasticity_param,&
								NTENS) 
		!Compute K_gamma_gamma
		A_gamma_gamma 	= -1.0_AbqRK
		K_gamma_gamma 	= A_gamma_Hc*B_Hc_gamma +A_gamma_gamma
		
		!Compute the X_gamma_epsilon quantity -> (-A_gamma_Hc*B_Hc_epsilon)/K_gamma_gamma
		X_gamma_epsilon = (-1.0_AbqRK*A_gamma_Hc*B_Hc_epsilon)/K_gamma_gamma
		
		!Compute partial derivate of stress by gamma (B_sigam_gamma)
		B_sigma_gamma 	= compute_derivative_stress_gamma(gamma_c_updated,C_0,T_c,stress_updated,&
														Strain_total_n2,Strain_inelastic_updated,&
														Inelasticity_param,NTENS)
		!Compute the material tangent
		C_t = B_sigma_epsilon + matmul(reshape(B_sigma_gamma,(/NTENS,1/))&
										,reshape(X_gamma_epsilon,(/1,NTENS/)))
		
	End Subroutine
	
	
	
	!--------------linearization functions and derivatives----------------------
	
	Function compute_A0_gamma(gamma_c,T_c,NTENS)result(A0_gamma)
		!=======================================================================
		!Function to compute strain localisation tensor A_0(gamma)
		!
		!	Inputs:-
		!	gama_c		-> micro crack density
		!	T_c 		-> Interaction tensor
		!	NTENS 		-> Stress components number
		!	
		!
		!	Outputs:-
		!   A0_gamma	-> strain localisation tensor A_0
		!=======================================================================
		IMPLICIT NONE
		!Dummy variable declarations
		integer(kind=AbqIK),intent(in)::NTENS
		real(kind=AbqRK),intent(in)::gamma_c
		real(kind=AbqRK),intent(in)::T_c(:,:)
		real(kind=AbqRK)::A0_gamma(NTENS,NTENS)
		
		!Local variables
		real(kind=AbqRK)::TempHolderA_0(NTENS,NTENS) !Temp holder for strain localisation tensor
		real(kind=AbqRK)::Identity_tensor(NTENS,NTENS)	!Fourth order identity tensor in voigt form
		!Variable initializations
		A0_gamma		= 0.0_AbqRK
		Identity_tensor	= 0.0_AbqRK
		TempHolderA_0 	= 0.0_AbqRK
		
		!--------Compute strain localisation tensor of virgin medium------------
		Identity_tensor	= I_4()
		TempHolderA_0	= Identity_tensor + gamma_c*(T_c - Identity_tensor)
		!A_0 is the inverse of TempHolderA_0
		call invert_matrix(TempHolderA_0,A0_gamma)
	End Function
	
	Function compute_A_gamma_Hc(gamma_c_infinity,H_c,WB_P)result(A_gamma_Hc)
		!=======================================================================
		!Function to compute partial derivative of phi_gamma with Hc
		!
		!	Inputs:-
		!	H_c 				-> yield criteria
		!	WB_P 				-> Weibul law parameters
		!	gamma_c_infinity 	-> micro-crack saturation parameter
		!	
		!
		!	Outputs:-
		!   A_gamma_Hc			-> partial derivative of phi_gamma with Hc
		!=======================================================================
		IMPLICIT NONE
		!Dummy variable declarations
		real(kind=AbqRK),intent(in)::gamma_c_infinity
		real(kind=AbqRK),intent(in)::H_c
		real(kind=AbqRK),intent(in)::WB_P(:)
		real(kind=AbqRK)::A_gamma_Hc,H_c_temp
	
		!Variable initializations
		A_gamma_Hc = 0.0_AbqRK
		H_c_temp = H_c
		!Condition for macually brackets
		if (H_c <= 1.0_AbqRK) H_c_temp = 1.0_AbqRK
		
		A_gamma_Hc = gamma_c_infinity*(WB_P(2)/WB_P(1))* &
						 exp(-((H_c_temp-1)/WB_P(1))**WB_P(2)) *&
						 ((H_c_temp-1)/WB_P(1))**(WB_P(2)-1)
	End Function
	
	Function compute_lambda_inelastic_strain(stress,Inelasticity_param,NTENS)result(Lambda)
		!=======================================================================
		!Function to compute lambda inelastic strain from flow rule (F_4:Stress)/H_s
		!
		!	Inputs:-
		!	stress					-> Stress state
		!	Inelasticity_param 		-> Inelastici parameters
		!	NTENS 					-> Stress components number
		!	
		!
		!	Outputs:-
		!   Lambda					-> (F_4:Stress)/H_s
		!=======================================================================
		IMPLICIT NONE
		!Dummy variable declarations
		integer(kind=AbqIK),intent(in)::NTENS
		real(kind=AbqRK),intent(in)::stress(:)
		real(kind=AbqRK),intent(in)::Inelasticity_param(:)
		real(kind=AbqRK)::Lambda(NTENS)
		
		!Local variables
		real(kind=AbqRK)::H_s 		!Inelastic strain yield value sqrt(Stress:F_4:Stress)
		!Variable initializations
		Lambda = 0.0_AbqRK
		H_s = 0.0_AbqRK
		
		H_s = sqrt((stress(2)*Inelasticity_param(1))**2+(stress(4)*Inelasticity_param(2))**2)
		Lambda(2) = (stress(2)*(Inelasticity_param(1)**2))/H_s
		Lambda(4) = (stress(4)*(Inelasticity_param(2)**2))/H_s
	End Function
	
	Function compute_derivative_A0_gamma(gamma_c,T_c,NTENS)result(dA0_gamma)
		!=======================================================================
		!Function to compute the partial derivative of A_0(gamma) by gamma
		!
		!	Inputs:-
		!	gama_c	-> micro crack density
		!	T_c 	-> Interaction tensor
		!	NTENS 	-> Stress components
		!	
		!
		!	Outputs:-
		!   dA0_gamma	-> partial derivative of A_0(gamma) by gamma
		!=======================================================================
		IMPLICIT NONE
		!Dummy variable declarations
		integer(kind=AbqIK),intent(in)::NTENS
		real(kind=AbqRK),intent(in)::gamma_c
		real(kind=AbqRK),intent(in)::T_c(:,:)
		real(kind=AbqRK)::dA0_gamma(NTENS,NTENS)
		
		!Local variable declarations
		real(kind=AbqRK)::TempHolderdA0_gamma(NTENS,NTENS)	!Temp variable to hold intermediate value
		real(kind=AbqRK)::A_0(NTENS,NTENS) 					!Strain localisation tensor
		real(kind=AbqRK)::Identity_tensor(NTENS,NTENS)		!Fourth order identity tensor in voigt form
		
		!Variable initializations
		dA0_gamma			= 0.0_AbqRK
		A_0					= 0.0_AbqRK
		TempHolderdA0_gamma = 0.0_AbqRK
		Identity_tensor 	= 0.0_AbqRK
		
		!Compute A_0(gamma)
		A_0			= compute_A0_gamma(gamma_c,T_c,NTENS)
		
		!Compute dA0_gamma
		Identity_tensor = I_4()
		TempHolderdA0_gamma = -1.0_AbqRK*matmul(A_0,(T_c-Identity_tensor))
		dA0_gamma 			= matmul(TempHolderdA0_gamma,A_0)
	
	End Function
	
	Function compute_B_Hc_gamma(H_c_i,C_0,T_c,Stress_i,Strain_total_n2,&
								Strain_inelastic_i,gamma_c_i,R,Inelasticity_param,&
								NTENS) result(B_Hc_gamma)
		!=======================================================================
		!Subroutine to compute the scalar value B_Hc_gamma for i th iteration of
		!cutting plane loop. Value comes from the linerization done for the
		!cutting plane algorithm. For ease of calculation the formula is
		!split into part_A,part_B and part_C. 
		!B_Hc_gamma = part_A:C_0:(part_B +part_C)
		!
		!	Inputs:-
		!	C_0						-> Initial material stiffness
		!	Strain_inelastic_i		-> Inelastic strain @ iteration i
		!	gamma_c_i				-> Micro crack density @ iteration i
		!	Stress_i				-> Trial stress @ iteration i
		!	R						-> Threshold material parameters
		!	gamma_c_infinity		-> Saturation micro crack density parameter
		!	Strain_total_n2 		-> Total strain @ time step n+1
		!	T_c						-> Interaction tensor
		!	Inelasticity_param 		-> Inelastici parameters
		!	NTENS					-> Number of stress components
		!	
		!
		!	Outputs:-
		!   B_Hc_gamma				-> Partial derivate of Hc with gamma
		!=======================================================================
		IMPLICIT NONE
		!Dummy variable declarations
		real(kind=AbqRK),intent(in)::H_c_i
		real(kind=AbqRK),intent(in)::C_0(:,:)
		real(kind=AbqRK),intent(in)::T_c(:,:)
		real(kind=AbqRK),intent(in)::Stress_i(:)
		real(kind=AbqRK),intent(in)::Strain_total_n2(:)
		real(kind=AbqRK),intent(in)::Strain_inelastic_i(:)
		real(kind=AbqRK),intent(in)::gamma_c_i
		real(kind=AbqRK),intent(in)::R(:)
		real(kind=AbqRK),intent(in)::Inelasticity_param(:)
		integer(kind=AbqIK),intent(in)::NTENS
		real(kind=AbqRK)::B_Hc_gamma
		
		!Local variable declarations
		real(kind=AbqRK)::part_A(NTENS) 				! Part A of B_Hc_gamma
		real(kind=AbqRK)::part_B(NTENS,1) 				! Part B of B_Hc_gamma
		real(kind=AbqRK)::part_C(NTENS,1) 				! Part C of B_Hc_gamma
		real(kind=AbqRK)::H_4(2) 						! Fourth order tesnor H from yield criteria( only two non zero components considered)
		real(kind=AbqRK)::local_stress(2) 				! Local stress state (S22 and S12 only taken)
		real(kind=AbqRK)::A_0(NTENS,NTENS) 				! Strain localisation tensor
		real(kind=AbqRK)::dA0_gamma(NTENS,NTENS) 		! Partial derivative of A_0(gamma) by gamma
		real(kind=AbqRK)::TempHolderB_Hc_gamma(1,NTENS) ! Temp variable for intermediate results
		real(kind=AbqRK)::Lambda_inelastic(NTENS) 		! (F_4:Stress)/H_s
		

		!Variable initializations
		part_A 					= 0.0_AbqRK
		part_B 					= 0.0_AbqRK
		part_C 					= 0.0_AbqRK
		H_4 					= 0.0_AbqRK
		local_stress 			= 0.0_AbqRK
		Lambda_inelastic 		= 0.0_AbqRK
		A_0 					= 0.0_AbqRK
		dA0_gamma 				= 0.0_AbqRK
		B_Hc_gamma 				= 0.0_AbqRK
		TempHolderB_Hc_gamma	= 0.0_AbqRK
		
		!Calculate part A -> (H_4:sigma_0)/H_c, as only two non zero components 
		! are there for the Fourth order tesnor H_4, only this values are computed.
		H_4(1) = 1.0_AbqRK/(R(1)**2)
		H_4(2) = 1.0_AbqRK/(R(2)**2)
		!Calculate the local stress components S22 and S12
		local_stress(1) = (1.0_AbqRK/(1.0_AbqRK-gamma_c_i))*Stress_i(2)
		local_stress(2) = (1.0_AbqRK/(1.0_AbqRK-gamma_c_i))*Stress_i(4)
		part_A(2) 		= (H_4(1)*local_stress(1))/H_c_i
		part_A(4) 		= (H_4(2)*local_stress(2))/H_c_i
		
		!Calculate part B -> (-A_0 : Lambda_inealstic)
		A_0					= compute_A0_gamma(gamma_c_i,T_c,NTENS)
		Lambda_inelastic	= compute_lambda_inelastic_strain(Stress_i,Inelasticity_param,NTENS)
		part_B 				= -1.0_AbqRK*matmul(A_0,reshape(Lambda_inelastic,(/NTENS,1/))) ! 6 x 6 * 6 x 1
		
		! Calculate part C -> -A_0:(T_c - I_4):A_0:(E-E_s)
		dA0_gamma	= compute_derivative_A0_gamma(gamma_c_i,T_c,NTENS)
		part_C		= matmul(dA0_gamma,reshape((Strain_total_n2-Strain_inelastic_i),(/NTENS,1/)))
		
		!Calculate B_Hc_gamma -> part_A : C_0 : (part_B +part_C)
		TempHolderB_Hc_gamma	= matmul( reshape(part_A,(/1,NTENS/)), C_0)
		
		!Mat multiplication of (1 x 6) * (6 x 1) equivalent to dot product of 2 vectors. 
		B_Hc_gamma				= dot_product(TempHolderB_Hc_gamma(1,:) ,(part_B(:,1) + part_C(:,1)))
		
	End Function
	
	Function compute_derivative_stress_strain(gamma_c,C_0,T_c,NTENS)result(B_sigma_epsilon)
		!=======================================================================
		!Function to compute the partial derivative of stress by strain
		!
		!	Inputs:-
		!	gama_c	-> micro crack density
		!	T_c 	-> Interaction tensor
		!	C_0 	-> Initial stiffness
		!	NTENS 	-> Stress components
		!	
		!
		!	Outputs:-
		!   B_sigma_epsilon	-> partial derivative of stress by strain
		!=======================================================================
		IMPLICIT NONE
		!Dummy variable declarations
		integer(kind=AbqIK),intent(in)::NTENS
		real(kind=AbqRK),intent(in)::gamma_c
		real(kind=AbqRK),intent(in)::C_0(:,:)
		real(kind=AbqRK),intent(in)::T_c(:,:)
		real(kind=AbqRK)::B_sigma_epsilon(NTENS,NTENS)
		
		!Local variable declarations
		real(kind=AbqRK)::D_stiffness_reduction(NTENS,NTENS)
	
		
		!Variable initializations
		D_stiffness_reduction = 0.0_AbqRK
		B_sigma_epsilon = 0.0_AbqRK
		
		!Compute stiffness reduction
		call compute_stiffness_reduction(gamma_c,D_stiffness_reduction,C_0,T_c,NTENS)
		
		!Compute B_sigma_epsilon
		B_sigma_epsilon = C_0 - D_stiffness_reduction
	End Function
	
	Function compute_derivative_Hc_strain(gamma_c,stress,C_0,T_c,R,NTENS)result(B_Hc_epsilon)
		!=======================================================================
		!Function to compute the partial derivative of Hc by strain
		!
		!	Inputs:-
		!	gama_c	-> micro crack density
		!	stress 	-> Stress
		!	T_c 	-> Interaction tensor
		!	C_0 	-> Initial stiffness
		!	NTENS 	-> Stress components
		!	
		!
		!	Outputs:-
		!   B_Hc_epsilon	-> partial derivative of Hc by strain
		!=======================================================================
		IMPLICIT NONE
		!Dummy variable declarations
		integer(kind=AbqIK),intent(in)::NTENS
		real(kind=AbqRK),intent(in)::gamma_c
		real(kind=AbqRK),intent(in)::stress(NTENS)
		real(kind=AbqRK),intent(in)::C_0(:,:)
		real(kind=AbqRK),intent(in)::T_c(:,:)
		real(kind=AbqRK),intent(in)::R(:)
		real(kind=AbqRK)::B_Hc_epsilon(NTENS)
		
		!Local variable declarations
		real(kind=AbqRK)::A_0(NTENS,NTENS) 							! Strain localization tensor
		real(kind=AbqRK)::local_stress(2) 							! Local stress ( s22,s12)
		real(kind=AbqRK)::H_c 										! Damage activation criteria
		real(kind=AbqRK)::H_4(2) 									! Fourth order tensor from damamge activation criteria. Only 2 non zero values taken
		real(kind=AbqRK)::double_contract_H4_local_stress(NTENS)	! (H_4:Sigma0)/H_c
		real(kind=AbqRK)::TempHolderB_Hc_epsilon(1,NTENS) 			! Temp variable for intermediate results
	
		
		!Variable initializations
		A_0								= 0.0_AbqRK
		B_Hc_epsilon					= 0.0_AbqRK
		local_stress					= 0.0_AbqRK
		H_c								= 0.0_AbqRK
		H_4 							= 0.0_AbqRK
		double_contract_H4_local_stress	= 0.0_AbqRK
		TempHolderB_Hc_epsilon 			= 0.0_AbqRK
		
		!Calculate the local stress components S22 and S12
		local_stress(1) = (1.0_AbqRK/(1.0_AbqRK-gamma_c))*stress(2)
		local_stress(2) = (1.0_AbqRK/(1.0_AbqRK-gamma_c))*stress(4)
		
		!Calculate Hc
		call compute_damage_activation_criteria(H_c,stress,gamma_c,R)
		
		!Caculate (H_4:local_stress)/Hc
		H_4(1) = 1.0_AbqRK/(R(1)**2)
		H_4(2) = 1.0_AbqRK/(R(2)**2)
		double_contract_H4_local_stress(2) = (H_4(1)*local_stress(1))/H_c
		double_contract_H4_local_stress(4) = (H_4(2)*local_stress(2))/H_c
		
		!Compute A_0(gamma)
		A_0 = compute_A0_gamma(gamma_c,T_c,NTENS)
		
		!Compute B_Hc_epsilon (Tensor in voigt notation) -> (H_4:local_stress)/Hc : C_0 : A_0
		TempHolderB_Hc_epsilon	= matmul( reshape(double_contract_H4_local_stress,(/1,NTENS/)), C_0)
		B_Hc_epsilon			= reshape(matmul(TempHolderB_Hc_epsilon ,A_0),(/NTENS/))
	End Function
	
	Function compute_derivative_D_gamma(gamma_c,T_c,C_0,NTENS)result(dD_gamma)
		!=======================================================================
		!Function to compute the partial derivative of D(gamma) by gamma
		!
		!	Inputs:-
		!	gama_c		-> micro crack density
		!	T_c 		-> Interaction tensor
		!	C_0 		-> Initial stiffness
		!	NTENS 		-> Stress components
		!	
		!
		!	Outputs:-
		!   dD_gamma	-> partial derivative of D(gamma) by gamma
		!=======================================================================
		IMPLICIT NONE
		!Dummy variable declarations
		integer(kind=AbqIK),intent(in)::NTENS
		real(kind=AbqRK),intent(in)::gamma_c
		real(kind=AbqRK),intent(in)::T_c(:,:)
		real(kind=AbqRK),intent(in)::C_0(:,:)
		real(kind=AbqRK)::dD_gamma(NTENS,NTENS)
		
		!Local variable declarations
		real(kind=AbqRK)::TempHolderdD_gamma1(NTENS,NTENS)	! Temp variable for intermediate results
		real(kind=AbqRK)::TempHolderdD_gamma2(NTENS,NTENS)	! Temp variable for intermediate results
		real(kind=AbqRK)::A_0(NTENS,NTENS)					! Strain localisation tensor
		real(kind=AbqRK)::dA0_gamma(NTENS,NTENS)			! Partial derivative of A0 by gamma
		
		
		!Variable initializations
		dD_gamma 			= 0.0_AbqRK
		A_0 				= 0.0_AbqRK
		dA0_gamma 			= 0.0_AbqRK
		TempHolderdD_gamma1	= 0.0_AbqRK
		TempHolderdD_gamma2	= 0.0_AbqRK
		
		!Compute A_0(gamma)
		A_0 = compute_A0_gamma(gamma_c,T_c,NTENS)
		
		!Compute dA_0(gamma)
		dA0_gamma = compute_derivative_A0_gamma(gamma_c,T_c,NTENS)
		
		
		!Compute dD_gamma
		TempHolderdD_gamma1 = matmul(C_0,T_c)
		TempHolderdD_gamma2 = A_0 + gamma_c*dA0_gamma
		dD_gamma = matmul(TempHolderdD_gamma1,TempHolderdD_gamma2)
	End Function
	
	Function compute_derivative_stress_gamma(gamma_c,C_0,T_c,stress,&
									Strain_total,Strain_inelastic,&
									Inelasticity_param,NTENS)result(B_sigma_gamma)
		!=======================================================================
		!Function to compute the partial derivative of stress by gamma
		!
		!	Inputs:-
		!	gama_c				-> micro crack density
		!	T_c 				-> Interaction tensor
		!	C_0 				-> Initial stiffness
		!	stress 				-> Stress
		!	Strain_total 		-> Total strain
		!	Strain_inelastic 	-> Inelastic strain
		!	Inelasticity_param 	-> Inelasticity parameter
		!	NTENS 				-> Stress components
		!	
		!
		!	Outputs:-
		!   B_sigma_gamma		-> partial derivative of stress by gamma
		!=======================================================================
		IMPLICIT NONE
		!Dummy variable declarations
		integer(kind=AbqIK),intent(in)::NTENS
		real(kind=AbqRK),intent(in)::gamma_c
		real(kind=AbqRK),intent(in)::C_0(:,:)
		real(kind=AbqRK),intent(in)::T_c(:,:)
		real(kind=AbqRK),intent(in)::stress(NTENS)
		real(kind=AbqRK),intent(in)::Strain_total(:)
		real(kind=AbqRK),intent(in)::Strain_inelastic(:)
		real(kind=AbqRK),intent(in)::Inelasticity_param(:)
		real(kind=AbqRK)::B_sigma_gamma(NTENS)
		
		!Local variable declarations
		real(kind=AbqRK)::D_stiffness_reduction(NTENS,NTENS)	! Stiffness reduction tensor
		real(kind=AbqRK)::stiffness_reduction(NTENS,NTENS) 		! Net reduction in stiffness (C_0 -D)
		real(kind=AbqRK)::Lambda_inealstic(NTENS) 				! (F_4:Stress)/H_s
		real(kind=AbqRK)::TempHolderB_sigma_gamma1(NTENS,1) 	! Temp variable for intermediate results
		real(kind=AbqRK)::TempHolderB_sigma_gamma2(NTENS,2) 	! Temp variable for intermediate results
		real(kind=AbqRK)::dD_gamma(NTENS,NTENS)
	
		
		!Variable initializations
		D_stiffness_reduction 		= 0.0_AbqRK
		dD_gamma 					= 0.0_AbqRK
		stiffness_reduction			= 0.0_AbqRK
		B_sigma_gamma 				= 0.0_AbqRK
		Lambda_inealstic			= 0.0_AbqRK
		TempHolderB_sigma_gamma1	= 0.0_AbqRK
		TempHolderB_sigma_gamma2	= 0.0_AbqRK
		
		!Compute stiffness reduction tensor D(gamma)
		call compute_stiffness_reduction(gamma_c,D_stiffness_reduction,C_0,T_c,NTENS)
		
		!Compute overall reduction in stiffness C_0 - D(gamma)
		stiffness_reduction = (C_0 - D_stiffness_reduction)
		
		!Compute Lambda strain from flow rule of inelastic strain
		Lambda_inealstic= compute_lambda_inelastic_strain(stress,Inelasticity_param,NTENS)
		
		!Compute partial derivative of stiffness reduction tensor by gamma
		dD_gamma = compute_derivative_D_gamma(gamma_c,T_c,C_0,NTENS)
		
		!Compute B_sigma_gamma -> -[C_0 - D(gamma)]:Lambda_inealstic - dD_gamma:(Total_strain-Inelastic_strain)
		TempHolderB_sigma_gamma1	= matmul(-1.0_AbqRK*stiffness_reduction,reshape(Lambda_inealstic,(/NTENS,1/)))
		TempHolderB_sigma_gamma2	= matmul(dD_gamma,reshape((Strain_total-Strain_inelastic),(/NTENS,1/)))
		B_sigma_gamma				= reshape(TempHolderB_sigma_gamma1 - TempHolderB_sigma_gamma2,(/NTENS/))
		
	End Function
	
	
	!---------------------Mathematical functions--------------------------------
	
	Subroutine invert_matrix(Matrix,Matrix_inv)
		!=======================================================================
		!Subroutine to invert a (square)matrix(M x N). Uses routines DGETRF and
		!DGETRI from LAPACK. Computes LU decomposition of matrix first using DGETRF.
		!Uses the LU decomposition to find the inverse using DGETRI.
		!Ref : LAPACK documentation
		!
		!	Inputs:-
		!	M	-> Matrix to be inverted
		!
		!	Outputs:-
		!   M_inv	-> Inverted matrix
		!=======================================================================
		IMPLICIT NONE
		!Dummy variable declarations
		real(kind=AbqRK),intent(in)::Matrix(:,:)
		real(kind=AbqRK),intent(out)::Matrix_inv(size(Matrix,1),size(Matrix,2))
		
		!Local variables
		real(kind=AbqRK)::work(size(Matrix,1)) !LAPACK work array
		integer(kind=AbqIK)::M !Number of rows of matrix A
		integer(kind=AbqIK)::N ! Number of columns of matrix A
		integer(kind=AbqIK)::LWORK ! The dimension of the array WORK.
		integer(kind=AbqIK)::LDA ! Leading dimension of A ie max(1,M)
		integer(kind=AbqIK)::info ! info returned by LAPACK routines
		integer(kind=AbqIK)::ipiv(size(Matrix,1)) !pivot indices from DGETRF. Dimension - min(M,N)
		real(kind=AbqRK)::LU(size(Matrix,1),size(Matrix,2))
		
		!Initializations
		M = size(Matrix,1)
		N = size(Matrix,2)
		LDA = size(Matrix,1)
		LWORK = M
		
		!First compute LU decomposition of matrix using DGETRF. 
		LU = Matrix
		call DGETRF(M,N,LU,LDA,ipiv,info)
		if(info .ne. 0) write(*,*) ' Matrix is singular. Not invertable  : '
		
		!If successfull, invert matrix using DGETRI. Uses LU decomposition and ipiv caluclated by DGETRF. Returns LU as the inverted matrix.
		call DGETRI(M,LU,LDA,ipiv,work,LWORK,info)
		if(info .ne. 0) write(*,*) ' Matrix is not inverted  : '
		
		!inverse
		Matrix_inv = LU
		
	End Subroutine
	
	Function I_4() result(Identity_4)
		!=======================================================================
		!Function to return the 4th order Identity tensor as a (6 x 6) matrix
		!	Inputs:-
		!	None
		!
		!	Outputs:-
		!   I_4	-> 4th order Identity tensor as a (6 x 6) matrix
		!=======================================================================
		IMPLICIT NONE
		!Dummy variable declarations
		real(kind=AbqRK)::Identity_4(6,6)
		
		!Local variable declarations
		integer(kind=AbqIK)::i,j
		
		!Variable initializations
		Identity_4(6,6) = 0.0_AbqRK
		
		!Assemble identity tensor
		forall(i=1:3,j=1:3,i==j)Identity_4(i,j)= 1.0_AbqRK
		forall(i=4:6,j=4:6,i==j)Identity_4(i,j)= 1.0_AbqRK
	
	End Function
	
	

END MODULE Yarn_utilities
