!DEC$ FREEFORM
!===============================================================================
! Seperate utility module for Matrix UMAT. Conatins all helper subroutines and
! functions called in UMAT Matrix.
!===============================================================================
MODULE Matrix_utilities
	USE ABQINTERFACE
	USE Type_State_variables_matrix
    USE Type_material_prop_matrix
    IMPLICIT NONE
    Private
	Public :: compute_viscoplastic_material_tangent,&
			  compute_StrainVe_Total, 			 	compute_stress,&
			  compute_viscoelastic_residual, 	 	compute_derivative_K_vi_vj,&
			  perform_viscoleastic_correction, 	 	compute_viscoelastic_material_tangent,&
			  compute_yield_function, 			 	compute_VonMises_eq_stress,&
			  compute_plastic_multiplier_residual, 	compute_damage_variable_residual,&
			  compute_total_energy_density_Y, 		compute_viscoplastic_residual,&
			  compute_derivative_K_p_vi, 			compute_derivative_K_r_vi,&
			  compute_derivative_K_D_vi, 			compute_derivative_K_r_r,&
			  compute_derivative_K_p_r, 			compute_derivative_K_D_p,&
			  assemble_derivatives_K_matrix,		perform_viscoplastic_correction
			  
			  
	Contains
	!-----------------------Procedures--------------------------------------
	
	Subroutine assemble_visco_elastic_stiffness(C_vi,PROPS,NTENS)
		!=======================================================================
		!Subroutine to assembel elastic stiffness matrix as 6*6 matrix from 
		!PROPS
		!	Inputs:-
		!	PROPS	-> Young's modulus and poisson's ratio
		!              
		!	Outputs:-
		!   C_0	-> Assembled 6x6 elastic stiffness in voigt form
		!=======================================================================
		IMPLICIT NONE
		!Dummy variable declarations
		real(kind=AbqRK),intent(in)::PROPS(:)
		integer(kind=AbqIK),intent(in)::NTENS
		real(kind=AbqRK),intent(out)::C_vi(NTENS,NTENS)
		
		!Local variables
		real(kind=AbqRK)::Lambda,mu
		real(kind=AbqRK)::E_e,nu
		integer(kind=AbqIK)::i,j
		
		!Variable initialization
		C_vi = 0.0_AbqRK
		E_e = PROPS(1)
		nu  = PROPS(2)
		
		!Compute Lame's constants : Lambda and mu
		Lambda =  E_e*nu/((1-2*nu)*(1+nu))
		mu = E_e/(2*(1+nu))
		
		!First half of diagonals
		forall(i=1:3)C_vi(i,i)= Lambda +2*mu
		!Second half of diagonals
		forall(i=4:6)C_vi(i,i)= mu
		!Off diagonal terms of upper matrix
		forall(i=1:3,j=1:3,i/=j)C_vi(i,j)= Lambda
	
	End Subroutine
	
	Subroutine compute_stress(mat_properties,Strain_total,Strain_Ve,Strain_Vp,Damage_var,Stress,NTENS)
		!=======================================================================
		!Subroutine to compute stress
		!
		!	Inputs:-
		!	C_e				-> Elastic stifness
		!	Strain_total_n2 -> Total strain @ n+1
		!   Strain_Ve_n2 	-> Viscoleastic strain from kelvin branches @ n +1    
		!	Strain_Vp_n 	-> Viscoplastic strain @ n
		!	N 				-> Number of kelvin branches 	 
		!    
		!	Outputs:-
		!   Stress_new		-> New stress value
		!=======================================================================
		IMPLICIT NONE
		!Dummy variable declarations
		Type(material_prop_matrix),intent(in)::mat_properties
		real(kind=AbqRK),intent(in)::Strain_total(:)
		real(kind=AbqRK),intent(in)::Strain_Ve(:,:),Strain_Vp(:)
		real(kind=AbqRK),intent(in)::Damage_var
		integer(kind=AbqIK),intent(in)::NTENS
		real(kind = AbqRK),intent(out)::Stress(NTENS)
		
		!Local variable declerations
		real(kind=AbqRK)::Strain_ve_total(NTENS)
		real(kind=AbqRK)::Strain(NTENS)
		
		!Variable initialization
		Stress	 		= 0.0_AbqRK
		Strain_ve_total	= 0.0_AbqRK
		Strain			= 0.0_AbqRK
		
		!Compute total viscoelastic strain from all Kelvin branches
		Strain_ve_total = compute_StrainVe_Total(Strain_Ve,mat_properties%N_k,NTENS)
		
		!Compute stress - (1-D)C_e : (Strain_total_n2-Strain_ve_total-Strain_Vp_n)
		Strain = Strain_total - Strain_ve_total - Strain_Vp
		Stress = (1-Damage_var)*reshape(MATMUL(mat_properties%C_e ,&
									   reshape(Strain,(/NTENS,1/))),(/NTENS/))
									   	
	End Subroutine
	
	!-------------------Routines to calculate resdiuals-------------------------
	
	Subroutine compute_viscoelastic_residual(Strain_total_n2,state_var_trial,state_var_previous,mat_properties,Ph_ve_i,DTIME,K_index,NTENS)
		!=======================================================================
		!Subroutine to compute stress
		!
		!	Inputs:-
		!	DTIME					-> delta time
		!	mat_properties					-> Elastic stifness
		!	K_index 			-> Index of kelvin branch
		!	Strain_total_n2 	-> Total strain @ n+1
		!   Strain_Ve_n2 		-> Viscoleastic strain from kelvin branches @ n+1 
		!	Outputs:-
		!    Ph_ve_i			-> Residual viscoelastic branch K_index
		!=======================================================================
		IMPLICIT NONE
		!Dummy variable declarations
		real(kind=AbqRK),intent(in)::Strain_total_n2(:)
		Type(state_variables_matrix),intent(in)::state_var_previous
		Type(material_prop_matrix),intent(in)::mat_properties
		Type(state_variables_matrix),intent(in)::state_var_trial
		real(kind=AbqRK),intent(in)::DTIME
		integer(kind=AbqIK),intent(in)::NTENS, K_index
		real(kind = AbqRK),intent(out)::Ph_ve_i(NTENS)
		
		!Local variable declerations
		real(kind=AbqRK)::Strain_ve_total(NTENS) 			!Total strain from all visco elastic branches
		real(kind=AbqRK)::Strain(NTENS) 					!Total strain ( E - sum(E_visco_eleastic) - E_visco_plastic)
		real(kind=AbqRK)::Strain_ve_delta(NTENS) 			!Change in viscoelastic strain E_vi@n+1 - E_vi@n
		real(kind=AbqRK)::V_inv_C_e,V_inv_Cv				!Quantity V_4_inv : C_e , V_4_inv : C_v
		integer::i
		!Variable initialization
		Strain_ve_total			= 0.0_AbqRK
		Strain					= 0.0_AbqRK
		Strain_ve_delta 		= 0.0_AbqRK
		V_inv_Cv 				= 0.0_AbqRK
		V_inv_C_e 				= 0.0_AbqRK
		Ph_ve_i					= 0.0_AbqRK
		
		!Compute delta viscoelastic strain
		Strain_ve_delta = reshape(((state_var_trial%Strain_Ve(K_index,:)-state_var_previous%Strain_Ve(K_index,:))/DTIME),(/NTENS/))
		
		!Compute total viscoelastic strain
		Strain_ve_total = compute_StrainVe_Total(state_var_trial%Strain_Ve,mat_properties%N_k,NTENS)
		
		!Compute total strain
		Strain = Strain_total_n2 - Strain_ve_total - state_var_trial%Strain_Vp
		
		!Compute V_4_inv : C_e
		V_inv_C_e = mat_properties%E_e/mat_properties%ViscoElastic_prop(K_index,2)
		
		!Compute V_inv:Cv 
		V_inv_Cv	= mat_properties%ViscoElastic_prop(K_index,1)/mat_properties%ViscoElastic_prop(K_index,2)
		
		!Compute the residual
		Ph_ve_i = Strain_ve_delta - V_inv_C_e*Strain + V_inv_Cv*state_var_trial%Strain_Ve(K_index,:)
			
	End Subroutine
	
	Subroutine compute_plastic_multiplier_residual(Strain_total_n2,state_var_trial,state_var_previous,mat_properties,DTIME,Phi_r,NDI,NTENS)
		!=======================================================================
		!Subroutine to compute resdiual of palstic multiplier		
		!Inputs:-
		!Plastic_mult_n2 	-> Plastic multiplier(r) @ n+1
		!Plastic_mult_n 		-> Plastic multiplier(r) @ n
		!DTIME				-> delta time
		!C_e					-> Elastic stifness
		!
		!Outputs:-
		!Phi_r				-> Residual of plastic multiplier r
		!=======================================================================
		IMPLICIT NONE
		!Dummy variable declarations
		Type(state_variables_matrix),intent(in)::state_var_previous
		Type(material_prop_matrix),intent(in)::mat_properties
		Type(state_variables_matrix),intent(in)::state_var_trial
		real(kind=AbqRK),intent(in)::Strain_total_n2(:), DTIME
		integer(kind=AbqIK),intent(in)::NTENS,NDI
		real(kind = AbqRK),intent(out)::Phi_r
		
		!Local variable declerations
		real(kind=AbqRK)::Intermediate_quantity1	!Intermediate result
		real(kind=AbqRK)::Intermediate_quantity2	!Intermediate result
		real(kind=AbqRK)::Intermediate_quantity3	!Intermediate result
		real(kind=AbqRK)::yield_function_value	    !Intermediate result
		real(kind = AbqRK)::Phi_r_FB
		real(kind=AbqRK)::H_p,K_p,R_0,n_exp,m_exp 	!Viscoplastic properties (Ref documentation for details)
		
		
		!Variable initialization
		Phi_r					= 0.0_AbqRK
		Intermediate_quantity1 	= 0.0_AbqRK
		Intermediate_quantity2 	= 0.0_AbqRK
		Intermediate_quantity3 	= 0.0_AbqRK
		R_0 					= mat_properties%ViscoPlastic_prop(1)
		K_p 					= mat_properties%ViscoPlastic_prop(2)
		n_exp 					= mat_properties%ViscoPlastic_prop(3)
		H_p 					= mat_properties%ViscoPlastic_prop(4)
		m_exp 					= mat_properties%ViscoPlastic_prop(5)
		
		!Compute Phi_r
		yield_function_value = compute_yield_function_positive(state_var_trial,mat_properties,NTENS,NDI)
		if(yield_function_value <= 1e-6_AbqRK) then
			Phi_r = 0.0_AbqRK
		else
			Intermediate_quantity1 = (state_var_trial%Plastic_mult-state_var_previous%Plastic_mult)/DTIME
			Intermediate_quantity2 = (1.0_AbqRK)/(H_p**(1.0_AbqRK/m_exp))
			Intermediate_quantity3 = (yield_function_value)**(1.0_AbqRK/m_exp)
		
			Phi_r = Intermediate_quantity1 - (Intermediate_quantity2 * Intermediate_quantity3)
		
		end if
			
	End Subroutine
	
	
	Subroutine compute_damage_variable_residual(Strain_total_n2,state_var_trial,state_var_previous,mat_properties,Phi_D,NTENS)
		!=======================================================================
		!Subroutine to compute resdiual of palstic multiplier
		!
		!	Inputs:-
		!	Plastic_mult_n 			-> Plastic multiplier(r) @ n
		!	Strain_total_n2 		-> Total strain @ n+1
		!	C_e						-> Elastic stifness 
		!    
		!	Outputs:-
		!   Phi_D 					-> Damage residual
		!    
		!=======================================================================
		IMPLICIT NONE
		!Dummy variable declarations
		Type(state_variables_matrix),intent(in)::state_var_previous
		Type(material_prop_matrix),intent(in)::mat_properties
		Type(state_variables_matrix),intent(in)::state_var_trial
		real(kind=AbqRK),intent(in)::Strain_total_n2(:)
		integer(kind=AbqIK),intent(in)::NTENS
		real(kind = AbqRK),intent(out)::Phi_D
		
		!Local variable declerations
		real(kind=AbqRK)::Y_total !Total energy density
		real(kind=AbqRK):: S_p,Beta ! Damage parameters
		
		!Variable initialization
		Y_total		= 0.0_AbqRK
		Phi_D		= 0.0_AbqRK
		S_p 	 	= mat_properties%ViscoPlastic_prop(6)
		Beta 		= mat_properties%ViscoPlastic_prop(7)
		
		!Compute energy density
		Y_total = compute_total_energy_density_Y(Strain_total_n2,state_var_trial,mat_properties,NTENS)
		
		!Compute Phi_D
		Phi_D = state_var_trial%Damage_var**2 - state_var_trial%Damage_var + state_var_previous%Damage_var - &
							state_var_previous%Damage_var*state_var_trial%Damage_var + &
							(((Y_total/S_p)**Beta)*(state_var_trial%Plastic_mult-state_var_previous%Plastic_mult))
			
	End Subroutine
	
	Subroutine compute_viscoplastic_residual(state_var_trial,state_var_previous,Phi_p,NDI,NTENS)
		!=======================================================================
		!Subroutine to compute viscoplastic residual
		!
		!	Inputs:-
		!	Strain_Vp_n2 		-> Viscoplastic strain @ n+1
		!	Strain_Vp_n 		-> Viscoplastic strain @ n
		!	Strain_total_n2 	-> Total strain @ n+1
		!    
		!	Outputs:-
		!   Phi_p 				-> Viscoplastic residual
		!    
		!=======================================================================
		IMPLICIT NONE
		!Dummy variable declarations
		Type(state_variables_matrix),intent(in)::state_var_previous
		Type(state_variables_matrix),intent(in)::state_var_trial
		integer(kind=AbqIK),intent(in)::NTENS,NDI
		real(kind = AbqRK),intent(out)::Phi_p(NTENS)
		
		!Local variable declerations
		real(kind=AbqRK)::Stress_dev(NTENS) !Deviatoric part of stress
		real(kind=AbqRK)::Stress_eq			!Von Mises equivalent stress
		real(kind=AbqRK)::Identity_2(NTENS) !2nd order Identity tensor
		
		!Variable initialization
		Phi_p		= 0.0_AbqRK
		Stress_dev	= 0.0_AbqRK
		Stress_eq	= 0.0_AbqRK
		Identity_2  = [1.0_AbqRK,1.0_AbqRK,1.0_AbqRK,0.0_AbqRK,0.0_AbqRK,0.0_AbqRK]
		
		!Compute the deviatoric part
		Stress_dev = state_var_trial%Stress - (1.0_AbqRK/3.0_AbqRK)*sum(state_var_trial%Stress(1:3))*Identity_2
		
		!Compute equivalent stress
		Stress_eq = compute_VonMises_eq_stress(state_var_trial%Stress,NDI,NTENS)
		
		!Compute residual
		Phi_p = state_var_trial%Strain_Vp -state_var_previous% Strain_Vp - 1.5_AbqRK*&
				((state_var_trial%Plastic_mult - state_var_previous%Plastic_mult)/(Stress_eq*(1-state_var_trial%Damage_var)))*Stress_dev
		
			
	End Subroutine
	
	!-----------Routines to perform visocelastic and palstic corrections--------
	
	Subroutine perform_viscoleastic_correction(Strain_total_n2,state_var_previous,state_var_trial,mat_properties,N,NTENS,DTIME)
		!=======================================================================
		!Subroutine to compute stress
		!
		!	Inputs:-
		!	DTIME	-> delta time
		!	C_e	-> Elastic stifness
		!	ViscoElastic_prop 	-> Viscoelastic properties of kelvin branches
		!	Strain_total_n2 -> Total strain @ n+1
		!   Strain_Ve_n2 -> Viscoleastic strain from kelvin branches @ n+1 
		!   Strain_Ve_n -> Viscoleastic strain from kelvin branches @ n     
		!	Strain_Vp_n -> Viscoplastic strain @ n
		!	Nu 		 -> poisson's ratio
		!	N -> Number of kelvin branches 	 
		!    
		!	Outputs:-
		!   Strain_Ve_n2	-> Corrected viscoleastic strains @ n+1
		!=======================================================================
		IMPLICIT NONE
		!Dummy variable declarations
		real(kind=AbqRK),intent(in),dimension(:)::Strain_total_n2
		real(kind=AbqRK),intent(in)::DTIME
		integer(kind=AbqIK),intent(in)::NTENS,N
		Type(state_variables_matrix),intent(in)::state_var_previous
		Type(material_prop_matrix),intent(in)::mat_properties
		Type(state_variables_matrix),intent(inout)::state_var_trial
		
		
		!Local variable declerations
		real(kind=AbqRK),parameter::TOEL 			 = 1e-6_AbqRK 	!Tolerence for convergence check
		integer(kind=AbqIK),parameter::Max_iteration = 20_AbqIK 	!Max number of iterations to break the loop
		real(kind = AbqRK)::delta_Strain_ve(N*NTENS) 				!Column matrix with all delta values
		real(kind = AbqRK)::Phi_ve(N*NTENS) 						!Column matrix with all residual values
		real(kind = AbqRK)::K_dd(N*NTENS,N*NTENS) 					!Coefficient matrix - with visco elastic derivatives
		real(kind = AbqRK)::delta_Strain_ve_matrix(N,NTENS) 		!For NR loop
		real(kind = AbqRK)::Phi_ve_i(NTENS) 						!Residulas for one kelvin branch
		real(kind = AbqRK)::K_ve_i_j(NTENS,NTENS) 					!Derivatives for one kelvin branch
		real(kind = AbqRK)::Residual_norm							!Norm of the residuals
		real(kind=AbqRK)::K_dd_LU(N*NTENS,N*NTENS)
		integer(kind=AbqIK)::ipiv(N*NTENS)
		integer(kind=AbqIK)::info 									! info returned by LAPACK routines
		integer(kind=AbqIK)::k,i,j,start_index,end_index,start_index2,end_index2

		
		!Variable initialization
		K_dd_LU 			= 0.0_AbqRK
		Residual_norm 		= 0.0_AbqRK
		delta_Strain_ve 	= 0.0_AbqRK
		K_dd				= 0.0_AbqRK
		Phi_ve				= 0.0_AbqRK
		K_ve_i_j			= 0.0_AbqRK
		
		!Initial guess @ k =1
		state_var_trial%Strain_Ve = state_var_previous%Strain_Ve
		state_var_trial%Strain_Vp = state_var_previous%Strain_Vp
			
		!compute and assemble residuals for all Kelvin branches
		start_index = 1
		do i=1,N
			Phi_ve_i = 0.0_AbqRK
			call compute_viscoelastic_residual(Strain_total_n2,state_var_trial,state_var_previous,mat_properties,Phi_ve_i,DTIME,i,NTENS)
			Phi_ve(start_index:) = -1.0_AbqRK*Phi_ve_i(:)
			start_index = start_index+NTENS
		end do
			
		!Compute and assemble the coeffient matrix with derivatives correspodning to each Kelvin branch
		start_index 	= 1
		start_index2 	= 1
		do i=1, N
			do j=1,N
				K_ve_i_j = compute_derivative_K_vi_vj(mat_properties,i,j,DTIME, NTENS)
				K_dd(start_index:,start_index2:) = K_ve_i_j(:,:)
				start_index2 = start_index2+size(K_ve_i_j,2)
				
			end do
				start_index2 = 1
				start_index = start_index+size(K_ve_i_j,1)
		end do
		
		!Factorize the coefficent matrix(A), for system of linear eqn A*X = B
		call LU_factorization_matrix(K_dd,K_dd_LU,ipiv)
		
		!solve the equation A*X = B and find the delta_Strain_ve values
		call dgetrs( 'N', size(K_dd_LU,1),1, K_dd_LU, size(K_dd_LU,1), ipiv, Phi_ve, size(Phi_ve,1), info )	
		
		if (info .eq. 0) then
			delta_Strain_ve 	 = Phi_ve
		else
			write(6,*) ' Linear system solver failed to solve for equation A*X = B. From matrix UMAT'
			delta_Strain_ve = 0.0_AbqRK
		end if 
		
		!Transform the column vecotr into matrix
		start_index = 1
		do i=1,N
			delta_Strain_ve_matrix(i,:) = delta_Strain_ve(start_index:)
			start_index = start_index+NTENS
		end do
			
		!update the visco elastic strain values
		state_var_trial%Strain_Ve 		= state_var_trial%Strain_Ve + delta_Strain_ve_matrix
		
		!set all other stae variables from previous time step
		state_var_trial%Strain_Vp 		= state_var_previous%Strain_Vp
		state_var_trial%Plastic_mult	= state_var_previous%Plastic_mult
		state_var_trial%Damage_var 		= state_var_previous%Damage_var
		
	End Subroutine
	
	
	Subroutine perform_viscoplastic_correction(Strain_total_n2,state_var_previous,state_var_trial,mat_properties,yield_f,DTIME,PNEWDT,NDI,N,NTENS,num_NR_iterations)
	!=======================================================================
	!Subroutine to compute stress
	!
	!	Inputs:-
	!	Strain_total_n2 	-> Total strain @ n+1
	!   Strain_Ve_n2 		-> Viscoleastic strain from kelvin branches @ n+1
	!   Strain_Ve_n 		-> Viscoleastic strain from kelvin branches @ n
	!	C_e					-> Elastic stifness
	!	ViscoElastic_prop	-> Viscoelastic properties
	!	K_dd_LU 			-> LU Decomposition of Derivative matrix for converged solution. Req for computing material tangent.
	!	Nu 		 			-> poisson's ratio
	!	DTIME				-> delta time 
	!    
	!	Outputs:-
	!   Strain_Vp_n2	-> Corrected viscoleastic strains @ n+1
	!=======================================================================
	IMPLICIT NONE
	!Dummy variable declarations
	real(kind=AbqRK),intent(in)::Strain_total_n2(:)
	real(kind=AbqRK),intent(in)::DTIME,yield_f
	integer(kind=AbqIK),intent(in)::NDI,N,NTENS
	integer(kind=AbqIK),intent(out)::num_NR_iterations
	real(kind=AbqRK),intent(inout)::PNEWDT
	Type(state_variables_matrix),intent(in)::state_var_previous
	Type(material_prop_matrix),intent(in)::mat_properties
	Type(state_variables_matrix),intent(inout)::state_var_trial
	
	
	
	!Local variable declarations
	real(kind=AbqRK),parameter::TOEL 			 = 0.5e-4_AbqRK 		!Tolerence for convergence check
	integer(kind=AbqIK),parameter::Max_iteration = 20_AbqIK 			!Max number of iterations to break the loop
	real(kind=AbqRK),dimension(6),parameter::Identity_tensor_2 = [1.0_AbqRK,1.0_AbqRK,1.0_AbqRK,0.0_AbqRK,0.0_AbqRK,0.0_AbqRK]
	Type(state_variables_matrix)::state_var_NR
	
	real(kind=AbqRK):: delta_Strain_vp(NTENS) 							!Delta strain values
	real(kind=AbqRK):: delta_Plastic_mult,delta_Damage_var 				!Delta damamge and palstic multiplier
	real(kind = AbqRK)::delta_Strain_ve_matrix(N,NTENS) 				!For NR loop
	real(kind = AbqRK)::Phi_column_matrix(N*NTENS+NTENS+2) 				!Column matrix with all residual values
	real(kind = AbqRK)::delta_column_matrix(N*NTENS+NTENS+2) 			!Column matrix with all delta unknows
	real(kind = AbqRK)::Phi_ve_i(NTENS),Phi_p(NTENS),Phi_r,Phi_D		!Residuals
	real(kind = AbqRK)::Residual_norm_Ve ,Residual_norm_Vp				!Norm of the residuals
	real(kind = AbqRK)::Residual_norm_r,Residual_norm_D					!Norm of the residuals
	real(kind = AbqRK)::TOEL_Ve,TOEL_r,TOEL_d,TOEL_vp
	real(kind = AbqRK)::K_dd(N*NTENS+NTENS+2,N*NTENS+NTENS+2)			!Matrix withh all partial derivatives
	real(kind=AbqRK)::K_dd_LU(N*NTENS+NTENS+2 ,N*NTENS+NTENS+2)
	real(kind=AbqRK)::Stress_proj(NTENS)								!Projected stress
	integer(kind=AbqIK)::ipiv(N*NTENS+NTENS+2),info 
	real(kind=AbqRK)::work_ve(N),work_vp(1),dlange
	integer(kind=AbqIK)::k,i,j
	integer(kind=AbqIK)::start_index,end_index,start_index2,end_index2
	logical::NR_Stabiliser_Flag
	
	!Allocate local NR state variables
	call state_var_NR % allocation(NTENS,N)
	
	!Initial guess @ k =1
	state_var_NR%Strain_Ve 		= state_var_trial%Strain_Ve
	state_var_NR%Strain_Vp 		= state_var_trial%Strain_Vp
	state_var_NR%Plastic_mult	= state_var_trial%Plastic_mult
	state_var_NR%Damage_var 	= state_var_trial%Damage_var
	state_var_NR%Stress 		= state_var_trial%Stress
	num_NR_iterations 			= 0.0_AbqIK
	Stress_proj 				= 0.0_AbqRK
	NR_Stabiliser_Flag 			= .True.
	
	!NR loop - start
	102 do k=1, Max_iteration
	
		!reset values
		K_dd_LU 				= 0.0_AbqRK
		K_dd					= 0.0_AbqRK
		Phi_column_matrix		= 0.0_AbqRK
		Phi_ve_i 				= 0.0_AbqRK
		Phi_p 					= 0.0_AbqRK
		Phi_r 					= 0.0_AbqRK
		Phi_D 					= 0.0_AbqRK
		delta_column_matrix 	= 0.0_AbqRK
		delta_Strain_ve_matrix 	= 0.0_AbqRK
		delta_Strain_vp 		= 0.0_AbqRK
		delta_Plastic_mult 		= 0.0_AbqRK
		delta_Damage_var 		= 0.0_AbqRK
		
		!In case of poor covergence stabilise local NR by improving the initial guess, using projected stress. Ref : Computational Plasticity, Page: 509
		if(k>=7 .and. NR_Stabiliser_Flag ==.True.) then
			Stress_proj = compute_projected_stress(state_var_trial,state_var_previous,mat_properties,NDI,NTENS)
			state_var_NR%Strain_Ve 		= state_var_trial%Strain_Ve
			state_var_NR%Strain_Vp 		= state_var_trial%Strain_Vp
			state_var_NR%Plastic_mult	= state_var_trial%Plastic_mult
			state_var_NR%Damage_var 	= state_var_trial%Damage_var
			state_var_NR%Stress 		= Stress_proj
			NR_Stabiliser_Flag = .False.
			go to 102
		end if
		
		!--compute and assemble all residual to the column matrix residuals -> Phi_column_matrix-------
		
		!compute and assemble visco elastic residuals for all Kelvin branches
		start_index = 1
		end_index 	= NTENS
		Residual_norm_Ve			= 0.0_AbqRK
		do i=1,N
			call compute_viscoelastic_residual(Strain_total_n2,state_var_NR,state_var_previous,mat_properties,Phi_ve_i,DTIME,i,NTENS)
			Phi_column_matrix(start_index:end_index) = -1.0_AbqRK*Phi_ve_i(:)
			start_index = start_index+NTENS
			end_index = end_index+NTENS
			!Compute norm of all residuls
			Residual_norm_Ve = Residual_norm_Ve+ sqrt(dot_product(Phi_ve_i(:),Phi_ve_i(:)))
		end do
		
		!compute and assemble visco plastic residual
		Residual_norm_Vp			= 0.0_AbqRK
		call compute_viscoplastic_residual(state_var_NR,state_var_previous,Phi_p,NDI,NTENS)
		Phi_column_matrix(N*NTENS+1:)= -1.0_AbqRK* Phi_p(:)
		Residual_norm_Vp = Residual_norm_Vp+ sqrt(dot_product(Phi_p(:),Phi_p(:)))

		!Compute and add plastic mutiplier residual
		Residual_norm_r = 0.0_AbqRK
		call compute_plastic_multiplier_residual(Strain_total_n2,state_var_NR,state_var_previous,mat_properties,DTIME,Phi_r,NDI,NTENS)
		Phi_column_matrix(N*NTENS+NTENS+1) =  -1.0_AbqRK* Phi_r
		Residual_norm_r =sqrt(Phi_r*Phi_r)
		
		!Compute and add damage variable residual
		Residual_norm_D = 0.0_AbqRK
		call compute_damage_variable_residual(Strain_total_n2,state_var_NR,state_var_previous,mat_properties,Phi_D,NTENS)
		Phi_column_matrix(N*NTENS+NTENS+2) =  -1.0_AbqRK* Phi_D
		Residual_norm_D = sqrt(Phi_D*Phi_D)
		
		!get tolerances
		TOEL_Ve	= TOEL*dlange('I',size(state_var_NR%Strain_Ve,1),size(state_var_NR%Strain_Ve,2),state_var_NR%Strain_Ve,size(state_var_NR%Strain_Ve,1), work_ve)
		TOEL_vp = TOEL*dlange('I',1,NTENS,reshape(state_var_NR%Strain_Vp,(/1,NTENS/)),1, work_vp)
		TOEL_r	= TOEL*abs(state_var_NR%Plastic_mult)
		TOEL_d  = TOEL*abs(state_var_NR%Damage_var)
		
		if ( Residual_norm_Ve > TOEL_Ve .or.  Residual_norm_Vp > TOEL_vp.or. &
			Residual_norm_r > TOEL_r .or.  Residual_norm_D > TOEL_d ) then
			
			!Compute and assemble the coeffient matrix with derivatives 
			K_dd = assemble_derivatives_K_matrix(Strain_total_n2,state_var_NR,state_var_previous,mat_properties,DTIME,NDI,N,NTENS)
											
			!Factorize the coefficent matrix(A), for system of linear eqn A*X = B
			call LU_factorization_matrix(K_dd,K_dd_LU,ipiv)
			
			!solve the equation A*X = B and find the delta_Strain_ve values
			call dgetrs( 'N', size(K_dd_LU,1),1, K_dd_LU, size(K_dd_LU,1), ipiv, Phi_column_matrix, size(Phi_column_matrix,1), info )	
			
			if (info .eq. 0) then
				delta_column_matrix 	 = Phi_column_matrix
			else
				write(6,*) ' Linear system solver failed to solve for equation A*X = B. From matrix UMAT'
				delta_column_matrix = 0.0_AbqRK
			end if 
			
			!Transform the column vector viscoelastic strain into matrix
			start_index = 1
			end_index 	= NTENS
			do i=1,N
				delta_Strain_ve_matrix(i,:) = delta_column_matrix(start_index:end_index)
				start_index = start_index+NTENS
				end_index = end_index+NTENS
			end do
			
			!Retrieve other delta strain values from the solution column vector
			delta_Strain_vp 	= delta_column_matrix(N*NTENS+1:)
			delta_Plastic_mult 	= delta_column_matrix(N*NTENS+NTENS+1)
			delta_Damage_var 	= delta_column_matrix(N*NTENS+NTENS+2)
			
			!update the visco elastic and plastic strain values and other state variables
			state_var_NR%Strain_Ve		= state_var_NR%Strain_Ve +  delta_Strain_ve_matrix
			state_var_NR%Strain_Vp		= state_var_NR%Strain_Vp+ delta_Strain_vp
			state_var_NR%Plastic_mult 	= state_var_NR%Plastic_mult+ delta_Plastic_mult
			state_var_NR%Damage_var		= state_var_NR%Damage_var + delta_Damage_var
			
			!Update the stress
			call compute_stress(mat_properties,Strain_total_n2,state_var_NR%Strain_Ve,state_var_NR%Strain_Vp,state_var_NR%Damage_var,state_var_NR%Stress,NTENS)
				
		else
			!update the trial state with the converged NR solutions
			state_var_trial%Strain_Ve 		= state_var_NR%Strain_Ve
			state_var_trial%Strain_Vp 		= state_var_NR%Strain_Vp
			state_var_trial%Plastic_mult	= state_var_NR%Plastic_mult
			state_var_trial%Damage_var 		= state_var_NR%Damage_var
			state_var_trial%Stress 			= state_var_NR%Stress
			num_NR_iterations 				= k
			if (k <=Max_iteration*0.7_AbqRK) PNEWDT=1.5_AbqRK
			EXIT
		
		end if
		
		if (k == Max_iteration ) then
			!If not converged, reset the state variables to previous state and request new time step.
			!write (6,*), 'Viscolplastic correction not converged.'
			PNEWDT=0.75_AbqRK
			state_var_trial%Strain_Ve 		= state_var_previous%Strain_Ve
			state_var_trial%Strain_Vp 		= state_var_previous%Strain_Vp
			state_var_trial%Plastic_mult	= state_var_previous%Plastic_mult
			state_var_trial%Damage_var 		= state_var_previous%Damage_var
			state_var_trial%Stress 			= state_var_previous%Stress
		end if
		
	end do
	
	!deallocate local NR state variables
	call state_var_NR % deallocation()
		
	End Subroutine
	
	!--------------routines to compute material tangents-------------------------
	
	Subroutine compute_viscoelastic_material_tangent(mat_properties,state_var_current,C_t,DTIME,N,NTENS)
		!=======================================================================
		!Subroutine to compute stress
		!
		!	Inputs:-
		!	C_e	-> Elastic stifness
		!	Damage_var_n		-> Damage variable @ n
		!
		!	Outputs:-
		!   C_t	-> Material tangent, when algorithm stops @ viscoelastic correction( ie no yielding)
		!=======================================================================
		IMPLICIT NONE
		!Dummy variable declarations
		Type(state_variables_matrix),intent(in)::state_var_current
		Type(material_prop_matrix),intent(in)::mat_properties
		real(kind=AbqRK),intent(in)::DTIME
		real(kind = AbqRK),intent(out)::C_t(NTENS,NTENS)
		integer(kind=AbqIK),intent(in)::NTENS,N
			
		
		!Local variable declerations
		real(kind = AbqRK)::K_ii(NTENS,NTENS) 					!Fourth order tensor with partial derivative of stress with strain
		real(kind = AbqRK)::K_sigma_vi(NTENS,NTENS) 			!Fourth order tensor - derivative of sigma with visco elastic strain
		real(kind = AbqRK)::K_id(NTENS,NTENS*N) 				!Fourth order tensor - derivative of residual with strain
		real(kind = AbqRK)::K_di(NTENS*N,NTENS) 				!Fourth order tensor - derivative of residual with strain
		real(kind = AbqRK)::K_di_red(NTENS*N,NTENS) 			!Reduced K_di. Obtained by soling Kdd * K_di_red = K_di
		real(kind = AbqRK)::K_vi_p(NTENS,NTENS) 				!Fourth order tensor - derivative of residual of viscoelastic strain with plastic strain
		real(kind=AbqRK)::Identity_tensor_4(NTENS,NTENS) 		!Fourth order identity tensor
		integer(kind=AbqIK)::i,j,start_index,start_index2		!Running indices
		real(kind = AbqRK)::K_dd_LU(N*NTENS,N*NTENS)
		real(kind = AbqRK)::K_dd(N*NTENS,N*NTENS)
		real(kind = AbqRK)::K_ve_i_j(NTENS,NTENS) 				!Derivatives for one kelvin branch
		integer(kind=AbqIK)::ipiv(N*NTENS) 
		integer(kind=AbqIK)::info 								!info returned by LAPACK routines
		
		!Variable initialization
		K_ii		= 0.0_AbqRK
		K_sigma_vi	= 0.0_AbqRK
		K_id		= 0.0_AbqRK
		K_di		= 0.0_AbqRK
		K_di_red	= 0.0_AbqRK
		K_vi_p		= 0.0_AbqRK
		K_dd_LU		= 0.0_AbqRK
		K_dd		= 0.0_AbqRK
		ipiv		= 0.0_AbqRK
		C_t 		= 0.0_AbqRK
		
		!Fourth order Identity tensor
		Identity_tensor_4 = I_4()
		
		!Compute and assemble the coeffient matrix with derivatives correspodning to each Kelvin branch
		start_index 	= 1
		start_index2 	= 1
		do i=1, N
			do j=1,N
				K_ve_i_j = compute_derivative_K_vi_vj(mat_properties,i,j,DTIME, NTENS)
				K_dd(start_index:,start_index2:) = K_ve_i_j(:,:)
				start_index2 = start_index2+size(K_ve_i_j,2)
				
			end do
				start_index2 = 1
				start_index = start_index+size(K_ve_i_j,1)
		end do
		
		!Factorize the coefficent matrix(A), for system of linear eqn A*X = B
		call LU_factorization_matrix(K_dd,K_dd_LU,ipiv)
		
		!Compute K_ii
		K_ii = (1-state_var_current%Damage_var)*mat_properties%C_e 
		
		!Assemble K_id
		start_index 	= 1
		start_index2 	= 1
		K_sigma_vi = -1.0_AbqRK*(1-state_var_current%Damage_var)*mat_properties%C_e
		do i=1, N
			K_id(start_index:,start_index2:) = K_sigma_vi(:,:)
			start_index2 = start_index2+size(K_sigma_vi,2)
		end do
		
		!Assemsble K_di
		start_index 	= 1
		do i=1, N
			K_vi_p = (mat_properties%E_e/mat_properties%ViscoElastic_prop(i,2))*Identity_tensor_4
			K_di(start_index:,:) = -1.0_AbqRK * K_vi_p(:,:)
			start_index = start_index+size(K_vi_p,1)
		end do
		
		
		!Solve and compute K_di_red = solve(K_dd_LU, RHS=K_di)
		call dgetrs( 'N', size(K_dd_LU,1),size(K_di,2), K_dd_LU, size(K_dd_LU,1), ipiv, K_di, size(K_di,1), info )	
		
		if (info .eq. 0) then
			K_di_red = K_di
		else
			write(6,*) ' Linear system solver failed to solve for equation A*X = B. From matrix UMAT-Material tangent viscoelastic correction',info
			K_di_red = 0.0_AbqRK
		end if 
		!Compute C_t
		C_t = K_ii - matmul(K_id,K_di_red)
		
		
	End Subroutine
	
	Subroutine compute_viscoplastic_material_tangent(Strain_total_n2,state_var_current,state_var_previous,mat_properties,C_t,NDI,DTIME,N,NTENS)
		!=======================================================================
		!Subroutine to compute material tangent when the algorithm enters 
		! viscoplastic corretion step
		!
		!	Inputs:-
		!	Stress_n2	-> Stress @ n+1
		!	Strain_total_n2	-> Total strain @ n+1
		!   Strain_Ve_n2 		-> Viscoleastic strain from kelvin branches @ (n+1)(k+1)
		!	Strain_Vp_n2 		-> Viscoplastic strain @ (n+1)(k+1)
		!	Plastic_mult_n2 	-> Plastic multiplier @ (n+1)(k+1)
		!	Plastic_mult_n	 	-> Plastic multiplier @ n
		!	C_e	-> Elastic stiffness
		!	ViscoElastic_prop	-> Viscoelastic properties
		!	ViscoPlastic_prop	-> Viscoelastic properties
		!	Nu	-> Poisson's ratio
		!    
		!	Outputs:-
		!   C_t	-> Material tangent, when algorithm stops @ viscoelastic correction( ie no yielding)
		!=======================================================================
		IMPLICIT NONE
		!Dummy variable declarations
		Type(state_variables_matrix),intent(in)::state_var_previous
		Type(material_prop_matrix),intent(in)::mat_properties
		Type(state_variables_matrix),intent(in)::state_var_current
		real(kind=AbqRK),intent(in)::Strain_total_n2(:)
		real(kind=AbqRK),intent(in)::DTIME
		integer(kind=AbqIK),intent(in)::NTENS,N,NDI
		real(kind = AbqRK),intent(out)::C_t(NTENS,NTENS)
		
		!Local variable declerations
		real(kind = AbqRK)::K_ii(NTENS,NTENS) 					! Deriavtive of stress with total strain
		real(kind = AbqRK)::K_di(N*NTENS+NTENS+2,NTENS) 		! Deriavtive of phi with total strain
		real(kind = AbqRK)::K_id(NTENS,N*NTENS+NTENS+2) 		! Deriavtive of Stress with state variables
		real(kind = AbqRK)::K_di_red(N*NTENS+NTENS+2,NTENS) 	! Reduced K_di. Obtained by soling Kdd * K_di_red = K_di
		real(kind = AbqRK)::K_dd(N*NTENS+NTENS+2,N*NTENS+NTENS+2)			!Matrix withh all partial derivatives
		real(kind=AbqRK)  ::K_dd_LU(N*NTENS+NTENS+2 ,N*NTENS+NTENS+2)
		integer(kind=AbqIK)::ipiv(N*NTENS+NTENS+2)
		integer(kind=AbqIK)::info,i ! info returned by LAPACK routines
		
		
		!Variable initialization
		C_t			= 0.0_AbqRK
		K_ii		= 0.0_AbqRK
		K_di		= 0.0_AbqRK
		K_id		= 0.0_AbqRK
		K_di_red	= 0.0_AbqRK
		ipiv		= 0.0_AbqRK
		K_dd		= 0.0_AbqRK
		
		!Compute and assemble the coeffient matrix with derivatives 
		K_dd = assemble_derivatives_K_matrix(Strain_total_n2,state_var_current,state_var_previous,mat_properties,DTIME,NDI,N,NTENS)
		
		!Factorize the coefficent matrix(A), for system of linear eqn A*X = B
		call LU_factorization_matrix(K_dd,K_dd_LU,ipiv)
		
		!Compute K_ii
		K_ii = (1-state_var_current%Damage_var)*mat_properties%C_e
		
		!Compute K_id
		K_id = compute_and_assemble_K_id(state_var_current,mat_properties,N,NTENS)
		
		
		!Compute K_di
		K_di = compute_and_assemble_K_di(Strain_total_n2,state_var_current,state_var_previous,mat_properties,N,NDI,DTIME,NTENS)
										
		!Solve and compute K_di_red = solve(K_dd_LU, RHS=K_di)
		!Call LAPACK subroutine DGESV to solve for X.
		call dgetrs( 'N', size(K_dd_LU,1),size(K_di,2), K_dd_LU, size(K_dd_LU,1), ipiv, K_di, size(K_di,1), info )	
		
		if (info .eq. 0) then
			K_di_red = K_di
		else
			write(6,*) ' Linear system solver failed to solve for equation A*X = B. From matrix UMAT-Material tangent viscoplastic correction'
			K_di_red = 0.0_AbqRK
		end if 
		
		!Compute C_t
		C_t = K_ii - matmul(K_id,K_di_red)
		
	End Subroutine
	
	
	!--------------functions to compute derivatives-----------------------------
	
	Function compute_StrainVe_Total(Strain_Ve_n,N,NTENS)result(Strain_ve_total)
		!=======================================================================
		!Function to compute total viscoelastic strain from N kelvin branches
		!
		!	Inputs:-
		!	Strain_Ve_n	-> Viscoelastic strain of Kelvin branches( each row corresponds to a branch)
		!	N	 		-> Number of Kelvin branches
		!	NTENS 		-> Stress components number
		!	
		!
		!	Outputs:-
		!   Strain_ve_total	-> Sum of viscoelastic strains from N branches
		!=======================================================================
		IMPLICIT NONE
		!Dummy variable declarations
		integer(kind=AbqIK),intent(in)::NTENS,N
		real(kind=AbqRK),intent(in)::Strain_Ve_n(:,:)
		real(kind=AbqRK)::Strain_ve_total(NTENS)
		
		!Local variables
		integer(kind=AbqIK)::i,j
		
		!Variable initializations
		Strain_ve_total		= 0.0_AbqRK
		
		!Compute the sum
		do i=1,N
			Strain_ve_total = Strain_ve_total + Strain_Ve_n(i,:)
		end do
	End Function
	
	Function  compute_derivative_K_vi_vj(mat_properties,Index_i,Index_j,DTIME, NTENS)result(K_vi_vj)
		!=======================================================================
		!Function to compute K_vi_vj
		!
		!	Inputs:-
		!	ViscoElastic_prop 	-> Viscoelastic properties of kelvin branches
		!	Outputs:-
		!   K_vi_vj				-> Derivative Ph_vi/Strain_ve_i
		!=======================================================================
		IMPLICIT NONE
		!Dummy variable declarations
		Type(material_prop_matrix),intent(in)::mat_properties
		real(kind=AbqRK),intent(in)::DTIME
		integer(kind=AbqIK),intent(in)::NTENS,Index_i,Index_j
		real(kind=AbqRK)::K_vi_vj(NTENS,NTENS)
		
		!Local variables
		real(kind=AbqRK)::Identity_tensor_4(NTENS,NTENS) !Fourth order identity tensor
		real(kind=AbqRK)::V_inv_Cv,V_inv_Ce		 		 !Quantity V^-1:C_vi
		
	
		!Variable initializations
		V_inv_Ce  			= 0.0_AbqRK
		V_inv_Cv 			= 0.0_AbqRK
		Identity_tensor_4 	= 0.0_AbqRK
		K_vi_vj 			= 0.0_AbqRK
		
		!Get the foruth order identity tensor
		Identity_tensor_4 = I_4()
		
		!Compute V_inv_Ce = V_vi_inv : C_e
		V_inv_Ce	= mat_properties%E_e/mat_properties%ViscoElastic_prop(Index_i,2)
		
		!Compute V_inv_Cv = V_vi_inv : C_vi
		V_inv_Cv	= mat_properties%ViscoElastic_prop(Index_i,1)/mat_properties%ViscoElastic_prop(Index_i,2)
		
		!Compute K_vi_vj
		if( Index_i == Index_j) then
		
			K_vi_vj = ((1.0_AbqRK/DTIME) + V_inv_Ce + V_inv_Cv) * Identity_tensor_4
		else
			K_vi_vj = V_inv_Ce*Identity_tensor_4
		
		end if
		
	End Function
	
	Function compute_derivative_K_vi_p(mat_properties,Index_i,NTENS)result(K_vi_p)
		!=======================================================================
		!Function to compute K_vi_p
		!
		!	Inputs:-
		!	ViscoElastic_prop 	-> Viscoelastic properties of kelvin branches
		!	Index_i,j 			-> Index of kelvin branch
		!	C_e					-> Elastic stifness
		!	Nu 		 			-> poisson's ratio
		!
		!	Outputs:-
		!   K_vi_p				-> Derivative Ph_vi/Strain_p
		!=======================================================================
		IMPLICIT NONE
		!Dummy variable declarations
		Type(material_prop_matrix),intent(in)::mat_properties
		integer(kind=AbqIK),intent(in)::NTENS,Index_i
		real(kind=AbqRK)::K_vi_p(NTENS,NTENS)
		
		!Local variables
		real(kind=AbqRK)::Identity_tensor_4(NTENS,NTENS) !Fourth order identity tensor
		real(kind=AbqRK)::V_4_inv_ce 					 !Quantity V_4_inv:C_e
		
	
		!Variable initializations
		V_4_inv_ce = 0.0_AbqRK
		K_vi_p = 0.0_AbqRK
		Identity_tensor_4 = 0.0_AbqRK
		
		!Get the foruth order identity tensor
		Identity_tensor_4 = I_4()
		
		!V_4_inv_ce
		V_4_inv_ce = mat_properties%E_e/mat_properties%ViscoElastic_prop(Index_i,2)
		
		!Compute K_vi_p
		K_vi_p	= V_4_inv_ce*Identity_tensor_4
		
	End Function
	
	Function compute_derivative_K_p_vi(state_var_trial,state_var_previous,mat_properties,NDI,NTENS) result(K_p_vi)
		!=======================================================================
		!Function to compute derivative K_p_vi(same for all i, so kelvin branch index not passed)
		!
		!	Inputs:-
		!	Plastic_mult_n2 	-> Plastic multiplier @ n+1
		!	Plastic_mult_n	 	-> Plastic multiplier @ n
		!	Stress_n2 			-> Stress @ n+1
		!	C_e					-> Elastic stifness
		!
		!	Outputs:-
		!   K_p_vi				-> Derivative Ph_p/Strain_ve_i
		!=======================================================================
		IMPLICIT NONE
		!Dummy variable declarations
		Type(state_variables_matrix),intent(in)::state_var_previous
		Type(material_prop_matrix),intent(in)::mat_properties
		Type(state_variables_matrix),intent(in)::state_var_trial
		integer(kind=AbqIK),intent(in)::NTENS,NDI
		real(kind=AbqRK)::K_p_vi(NTENS,NTENS)
		
		!Local variables
		real(kind=AbqRK),dimension(6),parameter::Identity_tensor_2 = [1.0_AbqRK,1.0_AbqRK,1.0_AbqRK,0.0_AbqRK,0.0_AbqRK,0.0_AbqRK]
		real(kind=AbqRK)::Identity_tensor_4(NTENS,NTENS) 		!Fourth order identity tensor
		real(kind=AbqRK)::P_dev_4(NTENS,NTENS) 					!Fourth order symm deviatoric projection operator
		real(kind=AbqRK)::C_e_dev(NTENS,NTENS) 					!Deviatoric part of C_e
		real(kind=AbqRK)::Stress_eq, Stress_dev(NTENS)  		!Von Mises equivalent stress and Deviatoric part of Stress
		real(kind=AbqRK)::Intermediate_quantity1 				!To store intermediate result
		real(kind=AbqRK)::Intermediate_quantity2(NTENS,NTENS) 	!To store intermediate result
		integer(kind=AbqIK)::i,j,k,l,m
		
		!Variable initializations
		Identity_tensor_4 		= 0.0_AbqRK
		C_e_dev 				= 0.0_AbqRK
		P_dev_4 				= 0.0_AbqRK
		K_p_vi 					= 0.0_AbqRK
		Intermediate_quantity1	= 0.0_AbqRK
		
		!Get the foruth order identity tensor
		Identity_tensor_4 = I_4()
		
		!Compute deviatoric projection P_dev_4
		P_dev_4 = Identity_tensor_4 - (1.0_AbqRK/3.0_AbqRK)*matmul(reshape(Identity_tensor_2,(/NTENS,1/)),reshape(Identity_tensor_2,(/1,NTENS/)))
		
		!Compute deviatoric part of stiffness
		!C_e_dev = 2*C_e(4,4)*P_dev_4
		C_e_dev = matmul(mat_properties%C_e,P_dev_4)
		
		!Compute the deviatoric part
		Stress_dev = state_var_trial%Stress - (1.0_AbqRK/3.0_AbqRK)*sum(state_var_trial%Stress(1:3))*Identity_tensor_2
		
		
		!Compute equivalent stress
		Stress_eq = compute_VonMises_eq_stress(state_var_trial%Stress,NDI,NTENS)
		
		!Compute -0.5*(r_n2-r_n)
		Intermediate_quantity1 = 1.5_AbqRK *(state_var_trial%Plastic_mult-state_var_previous%Plastic_mult)
		Intermediate_quantity2 = Intermediate_quantity1* ((Identity_tensor_4/Stress_eq) - &
								1.5_AbqRK*((matmul(reshape(Stress_dev,(/NTENS,1/)),reshape(Stress_dev,(/1,NTENS/))))/Stress_eq**3))
		
		K_p_vi = matmul(Intermediate_quantity2,C_e_dev)
		
	End Function
	
	Function compute_derivative_K_r_vi(state_var_trial,mat_properties,DTIME,NDI,NTENS)result(K_r_vi)
		!=======================================================================
		!Function to compute derivative K_r_vi(same for all i, so kelvin branch index not passed)
		!
		!	Inputs:-
		!	Stress_n2 			-> Stress @ n+1
		!	C_e					-> Elastic stifness
		!
		!	Outputs:-
		!   K_p_vi				-> Derivative Ph_p/Strain_ve_i
		!=======================================================================
		IMPLICIT NONE
		!Dummy variable declarations
		Type(material_prop_matrix),intent(in)::mat_properties
		Type(state_variables_matrix),intent(in)::state_var_trial
		real(kind=AbqRK),intent(in)::DTIME
		integer(kind=AbqIK),intent(in)::NTENS,NDI
		real(kind=AbqRK)::K_r_vi(NTENS)
		
		!Local variables
		real(kind=AbqRK),dimension(6),parameter::Identity_tensor_2 = [1.0_AbqRK,1.0_AbqRK,1.0_AbqRK,0.0_AbqRK,0.0_AbqRK,0.0_AbqRK]
		real(kind=AbqRK)::Identity_tensor_4(NTENS,NTENS) 		!Fourth order identity tensor
		real(kind=AbqRK)::P_dev_4(NTENS,NTENS) 					!Fourth order deviatoric projection operator
		real(kind=AbqRK)::C_e_dev(NTENS,NTENS) 					!Deviatoric part of C_e
		real(kind=AbqRK)::Stress_eq, Stress_dev(NTENS)			!Von Mises equivalent stress and Deviatoric part of Stress
		real(kind=AbqRK)::H_p,m_exp ,n_exp, K_p, R_0						!Viscoplastic properties (Ref documentation for details)
		real(kind=AbqRK)::Intermediate_quantity1 				!To store intermediate result
		real(kind=AbqRK)::Intermediate_quantity2 				!To store intermediate result
		real(kind=AbqRK)::yield_function_value 					!Positive value of yield function
		real(kind=AbqRK)::Intermediate_quantity3(NTENS) 		!To store intermediate result
		
		!Variable initializations
		Identity_tensor_4 		= 0.0_AbqRK
		C_e_dev 				= 0.0_AbqRK
		P_dev_4 				= 0.0_AbqRK
		K_r_vi 					= 0.0_AbqRK
		Intermediate_quantity2	= 0.0_AbqRK
		R_0						= mat_properties%ViscoPlastic_prop(1)
		K_p 					= mat_properties%ViscoPlastic_prop(2)
		n_exp 					= mat_properties%ViscoPlastic_prop(3)
		H_p 					= mat_properties%ViscoPlastic_prop(4)
		m_exp 					= mat_properties%ViscoPlastic_prop(5)
		
		
		!Get the foruth order identity tensor
		Identity_tensor_4 = I_4()
		
		!Compute deviatoric projection P_dev_4
		P_dev_4 = Identity_tensor_4 - (1.0_AbqRK/3.0_AbqRK)*matmul(reshape(Identity_tensor_2,(/NTENS,1/)),reshape(Identity_tensor_2,(/1,NTENS/)))
		
		!Compute deviatoric part of stiffness
		!C_e_dev =  2*C_e(4,4)*P_dev_4
		C_e_dev = matmul(mat_properties%C_e,P_dev_4)
		
		!Compute the deviatoric part
		Stress_dev = state_var_trial%Stress - (1.0_AbqRK/3.0_AbqRK)*sum(state_var_trial%Stress(1:3))*Identity_tensor_2
		
		
		!Compute equivalent stress
		Stress_eq = compute_VonMises_eq_stress(state_var_trial%Stress,NDI,NTENS)
		
		!Compute K_r_vi
		yield_function_value  = compute_yield_function_positive(state_var_trial,mat_properties,NTENS,NDI)
		if(yield_function_value <= 1e-6_AbqRK) then
			K_r_vi = 0.0_AbqRK
		else
			Intermediate_quantity1 = 1.5_AbqRK/(m_exp*(H_p**(1/m_exp))*Stress_eq)
			Intermediate_quantity2 = (yield_function_value)**((1-m_exp)/m_exp)
			Intermediate_quantity3 = Intermediate_quantity1 * Intermediate_quantity2 * Stress_dev
			K_r_vi = reshape(matmul(reshape(Intermediate_quantity3,(/1,NTENS/)),C_e_dev),(/NTENS/))
		end if
		
	End Function
	
	
	Function compute_derivative_K_D_vi(Strain_total_n2,state_var_trial,state_var_previous,mat_properties,K_index,N,NTENS)result(K_D_vi)
		!=======================================================================
		!Function to compute derivative K_D_vi(for i'th (K_index) kelvin branch)
		!(To avoid ambiguity, the large eqn is broken down into intermediate results)
		!
		!	Inputs:-
		!	Strain_total_n2 	-> Total strain @ n+1
		!	C_e					-> Elastic stifness
		!	ViscoElastic_prop	-> Viscoelastic properties
		!
		!	Outputs:-
		!   K_D_vi				-> Derivative Phi_D/E_vi
		!=======================================================================
		IMPLICIT NONE
		!Dummy variable declarations
		Type(state_variables_matrix),intent(in)::state_var_previous
		Type(material_prop_matrix),intent(in)::mat_properties
		Type(state_variables_matrix),intent(in)::state_var_trial
		real(kind=AbqRK),intent(in)::Strain_total_n2(:)
		integer(kind=AbqIK),intent(in)::NTENS,N,K_index
		real(kind=AbqRK)::K_D_vi(NTENS)
		
		!Local variables
		real(kind=AbqRK)::Intermediate_quantity1           !To store intermediate result
		real(kind=AbqRK)::Intermediate_quantity2(NTENS)    !To store intermediate result
		real(kind=AbqRK)::Strain_elastic_n2(NTENS)         !Total elastic strain @ n+1
		real(kind=AbqRK)::Intermediate_quantity3(NTENS)    !To store intermediate result
		real(kind=AbqRK)::Strain_Ve_total(NTENS) 		   !Total viscoelastic strain
		real(kind=AbqRK)::C_vi(NTENS,NTENS),Strain_Ve_n2_i !Stiffness of viscoelastic spring i, Viscoelastic strain of i'the Kelvin branch
		real(kind=AbqRK)::S_p,Beta 				 		   !Damage parameters
		real(kind=AbqRK)::Y_total				 		   !Total energy density
		
	
		!Variable initializations
		K_D_vi					= 0.0_AbqRK
		Y_total					= 0.0_AbqRK
		Intermediate_quantity1	= 0.0_AbqRK
		Strain_elastic_n2 		= 0.0_AbqRK
		Intermediate_quantity3	= 0.0_AbqRK
		Strain_Ve_total 		= 0.0_AbqRK
		S_p 					= mat_properties%ViscoPlastic_prop(6)
		Beta 					= mat_properties%ViscoPlastic_prop(7)
		
		!Compute Total energy density
		Y_total = compute_total_energy_density_Y(Strain_total_n2,state_var_trial,mat_properties,NTENS)
				                  
		!Compute intermediate result 1 -> beta*(Y/S)**(beta-1)*(Plastisc_mul_n2 - Plastic_mult_n)
		Intermediate_quantity1 = (Beta/S_p)*((Y_total/S_p)**(Beta-1.0_AbqRK))*(state_var_trial%Plastic_mult-state_var_previous%Plastic_mult)
		

		!Compute Intermediate result2 -> -Stress_n2/(1-D)
		Intermediate_quantity2 = (-1.0_AbqRK * state_var_trial%Stress)/(1-state_var_trial%Damage_var)
		
		!Compute Stiffness of viscoelastic spring i
		call assemble_visco_elastic_stiffness(C_vi,[mat_properties%ViscoElastic_prop(K_index,1),mat_properties%Nu],NTENS)
		
		!Compute Intermediate_quantity3 -> C_vi : Strain_Ve_ni
		Intermediate_quantity3 = reshape(matmul(C_vi,reshape(state_var_trial%Strain_Ve(K_index,:),(/NTENS,1/))) ,(/NTENS/))
		
		!Compute K_D_vi -> Intermediate_quantity1 * (Intermediate_quantity2 + Intermediate_quantity3)
		K_D_vi = Intermediate_quantity1 * (Intermediate_quantity2 + Intermediate_quantity3)
		
	End Function
	
	Function  compute_derivative_K_r_r(state_var_trial,mat_properties,DTIME,NDI,NTENS)result(K_r_r)
		!=======================================================================
		!Function to compute derivative K_r_r
		!
		!	Inputs:-
		!	Plastic_mult_n2 	-> Plastic multiplier @ n+1
		!	ViscoPlastic_prop	-> Viscoplastic properties
		!
		!	Outputs:-
		!   K_r_r				-> Derivative Ph_r/r
		!=======================================================================
		IMPLICIT NONE
		!Dummy variable declarations
		Type(material_prop_matrix),intent(in)::mat_properties
		Type(state_variables_matrix),intent(in)::state_var_trial
		real(kind=AbqRK),intent(in)::DTIME
		integer(kind=AbqIK),intent(in)::NTENS,NDI
		real(kind=AbqRK)::K_r_r
		
		!Local variables
		real(kind=AbqRK)::Stress_eq				 	!Von Mises equivalent stress
		real(kind=AbqRK)::Intermediate_quantity1 	!To store intermediate result
		real(kind=AbqRK)::Intermediate_quantity2 	!To store intermediate result
		real(kind=AbqRK)::Intermediate_quantity3 	!To store intermediate result
		real(kind=AbqRK)::yield_function_value 		!To store intermediate result
		real(kind=AbqRK)::R_0,K_p,H_p,m_exp, n_exp	!Viscoplastic properties (Ref documentation for details)
		
	
		!Variable initializations
		K_r_r = 0.0_AbqRK
		Intermediate_quantity1	= 0.0_AbqRK
		Intermediate_quantity2	= 0.0_AbqRK
		Intermediate_quantity3	= 0.0_AbqRK
		R_0 					= mat_properties%ViscoPlastic_prop(1)
		K_p 					= mat_properties%ViscoPlastic_prop(2)
		n_exp					= mat_properties%ViscoPlastic_prop(3)
		H_p 					= mat_properties%ViscoPlastic_prop(4)
		m_exp 					= mat_properties%ViscoPlastic_prop(5)
		
		!Compute equivalent stress
		Stress_eq = compute_VonMises_eq_stress(state_var_trial%Stress,NDI,NTENS)
		
		!Compute K_r_r : split into 3 intermediate results
		yield_function_value   = compute_yield_function_positive(state_var_trial,mat_properties,NTENS,NDI)
		
		if(yield_function_value <=1e-6_AbqRK ) then
			K_r_r = (1.0_AbqRK/DTIME)
		else
		
			Intermediate_quantity1 = n_exp/(m_exp*H_p**(1/m_exp))
			Intermediate_quantity2 = (yield_function_value)**((1-m_exp)/m_exp)
			
			if(state_var_trial%Plastic_mult <= 1e-6_AbqRK) then
				Intermediate_quantity3 = K_p*1e-6_AbqRK**(n_exp-1)
			else
				Intermediate_quantity3 = K_p*state_var_trial%Plastic_mult**(n_exp-1)
			end if
			
			K_r_r = (1.0_AbqRK/DTIME)+ (Intermediate_quantity1*Intermediate_quantity2 * Intermediate_quantity3)
		end if
		
	End Function
	
	Function compute_derivative_K_p_D(state_var_trial,state_var_previous,mat_properties,DTIME,NDI,NTENS)result(K_p_D)
		!=======================================================================
		!Function to compute derivative K_p_D
		!
		!	Inputs:-
		!	Plastic_mult_n2 	-> Plastic multiplier @ n+1
		!	ViscoPlastic_prop	-> Viscoplastic properties
		!
		!	Outputs:-
		!   K_p_D				-> Derivative Ph_p/D
		!=======================================================================
		IMPLICIT NONE
		!Dummy variable declarations
		Type(state_variables_matrix),intent(in)::state_var_previous
		Type(material_prop_matrix),intent(in)::mat_properties
		Type(state_variables_matrix),intent(in)::state_var_trial
		real(kind=AbqRK),intent(in)::DTIME
		integer(kind=AbqIK),intent(in)::NTENS,NDI
		real(kind=AbqRK)::K_p_D(NTENS)
		
		!Local variables
		real(kind=AbqRK),dimension(6),parameter::Identity_tensor_2 = [1.0_AbqRK,1.0_AbqRK,1.0_AbqRK,0.0_AbqRK,0.0_AbqRK,0.0_AbqRK]
		real(kind=AbqRK)::Identity_tensor_4(NTENS,NTENS) 		!Fourth order identity tensor
		real(kind=AbqRK)::P_dev_4(NTENS,NTENS) 					!Fourth order deviatoric projection operator
		real(kind=AbqRK)::C_e_dev(NTENS,NTENS) 					!Deviatoric part of C_e
		real(kind=AbqRK)::Stress_eq, Stress_dev(NTENS)			!Von Mises equivalent stress and Deviatoric part of Stress
		real(kind=AbqRK)::Intermediate_quantity1 	!To store intermediate result
		real(kind=AbqRK)::yield_function_value 	!To store intermediate result
		real(kind=AbqRK)::R_0,K_p,H_p,m_exp, n_exp	!Viscoplastic properties (Ref documentation for details)
		
	
		!Variable initializations
		K_p_D = 0.0_AbqRK
		Intermediate_quantity1	= 0.0_AbqRK
		R_0 					= mat_properties%ViscoPlastic_prop(1)
		K_p 					= mat_properties%ViscoPlastic_prop(2)
		n_exp					= mat_properties%ViscoPlastic_prop(3)
		H_p 					= mat_properties%ViscoPlastic_prop(4)
		m_exp 					= mat_properties%ViscoPlastic_prop(5)
		
		!Compute the deviatoric part
		Stress_dev = state_var_trial%Stress - (1.0_AbqRK/3.0_AbqRK)*sum(state_var_trial%Stress(1:3))*Identity_tensor_2
		
		!Compute equivalent stress
		Stress_eq = compute_VonMises_eq_stress(state_var_trial%Stress,NDI,NTENS)
		
		Intermediate_quantity1 = -1.5_AbqRK *((state_var_trial%Plastic_mult-state_var_previous%Plastic_mult)/(Stress_eq*(1-state_var_trial%Damage_var)**2))
		K_p_D = Intermediate_quantity1 * (Stress_dev)
		
	End Function
	
	Function compute_derivative_K_p_r(state_var_trial,NDI,NTENS)result(K_p_r)
		!=======================================================================
		!Function to compute derivative K_p_r
		!
		!	Inputs:-
		!	Stress_n2 		-> Stress @ n+1
		!   Damage_var_n2 	-> Damage variable @ n+1      
		!
		!	Outputs:-
		!   K_p_r			-> Derivative Phi_p/r
		!=======================================================================
		IMPLICIT NONE
		!Dummy variable declarations
		Type(state_variables_matrix),intent(in)::state_var_trial
		integer(kind=AbqIK),intent(in)::NTENS,NDI
		real(kind=AbqRK)::K_p_r(NTENS)
		
		!Local variables
		real(kind=AbqRK),dimension(6),parameter::Identity_tensor_2 = [1.0_AbqRK,1.0_AbqRK,1.0_AbqRK,0.0_AbqRK,0.0_AbqRK,0.0_AbqRK]
		real(kind=AbqRK)::Stress_dev(NTENS),Stress_eq !Deviatoric part of stress and VM Stress
		
		!Variable initializations
		K_p_r = 0.0_AbqRK
		Stress_dev 			= 0.0_AbqRK
		Stress_eq 			= 0.0_AbqRK
		
		!Compute the deviatoric part
		Stress_dev = state_var_trial%Stress - (1.0_AbqRK/3.0_AbqRK)*sum(state_var_trial%Stress(1:3))*Identity_tensor_2
		
		
		!Compute equivalent stress
		Stress_eq = compute_VonMises_eq_stress(state_var_trial%Stress,NDI,NTENS)
		
		!Compute K_p_r
		K_p_r = (-1.5_AbqRK/(Stress_eq*(1-state_var_trial%Damage_var)))*Stress_dev
		
	End Function
	
	Function compute_derivative_K_D_p(Strain_total_n2,state_var_trial,state_var_previous,mat_properties,N,NTENS)result(K_D_p)
		!=======================================================================
		!Function to compute derivative K_D_vi(for i'th (K_index) kelvin branch)
		!(To avoid ambiguity, the large eqn is broken down into intermediate results)
		!
		!	Inputs:-
		!	Strain_total_n2 	-> Total strain @ n+1
		!	ViscoPlastic_prop	-> Viscoplastic properties
		!
		!	Outputs:-
		!   K_D_p				-> Derivative Phi_D/E_p
		!=======================================================================
		IMPLICIT NONE
		!Dummy variable declarations
		Type(state_variables_matrix),intent(in)::state_var_previous
		Type(material_prop_matrix),intent(in)::mat_properties
		Type(state_variables_matrix),intent(in)::state_var_trial
		real(kind=AbqRK),intent(in)::Strain_total_n2(:)
		integer(kind=AbqIK),intent(in)::NTENS,N
		real(kind=AbqRK)::K_D_p(NTENS)
		
		!Local variables
		real(kind=AbqRK)::Intermediate_quantity1 	!To store intermediate result
		real(kind=AbqRK)::Strain_Ve_total(NTENS) 	!To store intermediate result
		real(kind=AbqRK)::Strain_elastic_n2(NTENS)  !To store intermediate result
		real(kind=AbqRK)::S_p,Beta 				 	!Damage parameters
		real(kind=AbqRK)::Y_total				 	!Total energy density
		
	
		!Variable initializations
		K_D_p 					= 0.0_AbqRK
		Y_total 				= 0.0_AbqRK
		Strain_Ve_total 		= 0.0_AbqRK
		Strain_elastic_n2 		= 0.0_AbqRK
		Intermediate_quantity1 	= 0.0_AbqRK
		S_p 					= mat_properties%ViscoPlastic_prop(6)
		Beta 					= mat_properties%ViscoPlastic_prop(7)
		
		!Compute Total energy density
		Y_total = compute_total_energy_density_Y(Strain_total_n2,state_var_trial,mat_properties,NTENS)
				                  
		!Compute intermediate result 1 -> beta*(Y/S)**(beta-1)*(Plastisc_mul_n2 - Plastic_mult_n)
		Intermediate_quantity1 = -1.0_AbqRK*(Beta/S_p)*((Y_total/S_p)**(Beta-1.0_AbqRK))*(state_var_trial%Plastic_mult-state_var_previous%Plastic_mult)
		
		!Compute  K_D_p
		K_D_p = (Intermediate_quantity1/(1-state_var_trial%Damage_var))*state_var_trial%Stress
		
		
	End Function
	
	Function compute_derivative_K_r_D(state_var_trial,state_var_previous,mat_properties,DTIME,NDI,NTENS)result(K_r_D)
		!=======================================================================
		!Function to compute derivative K_r_D
		!
		!	Inputs:-
		!	Plastic_mult_n2 	-> Plastic multiplier @ n+1
		!	ViscoPlastic_prop	-> Viscoplastic properties
		!
		!	Outputs:-
		!   K_r_D				-> Derivative Ph_r/D
		!=======================================================================
		IMPLICIT NONE
		!Dummy variable declarations
		Type(state_variables_matrix),intent(in)::state_var_previous
		Type(material_prop_matrix),intent(in)::mat_properties
		Type(state_variables_matrix),intent(in)::state_var_trial
		real(kind=AbqRK),intent(in)::DTIME
		integer(kind=AbqIK),intent(in)::NTENS,NDI
		real(kind=AbqRK)::K_r_D
		
		!Local variables
		real(kind=AbqRK)::Stress_eq				 	! Von Mises equivalent stress
		real(kind=AbqRK)::Intermediate_quantity1 	!To store intermediate result
		real(kind=AbqRK)::Intermediate_quantity2 	!To store intermediate result
		real(kind=AbqRK)::Intermediate_quantity3 	!To store intermediate result
		real(kind=AbqRK)::yield_function_value 	!To store intermediate result
		real(kind=AbqRK)::R_0,K_p,H_p,m_exp, n_exp	!Viscoplastic properties (Ref documentation for details)
		
	
		!Variable initializations
		K_r_D = 0.0_AbqRK
		Intermediate_quantity1	= 0.0_AbqRK
		Intermediate_quantity2	= 0.0_AbqRK
		Intermediate_quantity3	= 0.0_AbqRK
		R_0 					= mat_properties%ViscoPlastic_prop(1)
		K_p 					= mat_properties%ViscoPlastic_prop(2)
		n_exp					= mat_properties%ViscoPlastic_prop(3)
		H_p 					= mat_properties%ViscoPlastic_prop(4)
		m_exp 					= mat_properties%ViscoPlastic_prop(5)
		
		!Compute equivalent stress
		Stress_eq = compute_VonMises_eq_stress(state_var_trial%Stress,NDI,NTENS)
		
		!Compute K_r_r : split into 3 intermediate results
		yield_function_value   = compute_yield_function_positive(state_var_trial,mat_properties,NTENS,NDI)
		
		if(yield_function_value <= 1e-6_AbqRK) then
			K_r_D  = 0.0_AbqRK
		else
			Intermediate_quantity1 = (1.0_AbqRK/(H_p**(1.0_AbqRK/m_exp)))
			Intermediate_quantity2 = (yield_function_value)**((1-m_exp)/m_exp)
			Intermediate_quantity3 = (Stress_eq)/(1.0_AbqRK-state_var_trial%Damage_var)**2
			
			K_r_D = Intermediate_quantity1 *Intermediate_quantity2*Intermediate_quantity3
		end if
		
	End Function
	
	
	!--------------------------Misc functions-----------------------------------
	
			 
	Function compute_yield_function(state_var_trial,state_var_previous,mat_properties,NTENS,NDI)result(yield_f)
		!=======================================================================
		!Function to compute the value of yield function (f)
		!
		!	Inputs:-
		!	Stress_n2			-> updated stress computed with the converged viscoelastic strain value
		!	ViscoPlastic_prop	-> Viscoelastic properties of kelvin branches
		!	NTENS				-> Stress components number
		!	Damage_var_n		-> Damage variable @ n
		!	Plastic_mult_n		-> Plastic multiplier @ n
		!	
		!
		!	Outputs:-
		!   yield_f				-> yield function value
		!=======================================================================
		IMPLICIT NONE
		!Dummy variable declarations
		Type(state_variables_matrix),intent(in)::state_var_previous
		Type(material_prop_matrix),intent(in)::mat_properties
		Type(state_variables_matrix),intent(in)::state_var_trial
		integer(kind=AbqIK),intent(in)::NTENS,NDI
		real(kind=AbqRK)::yield_f
		
		!Local variables
		integer(kind=AbqIK)::i,j
		real(kind=AbqRK)::K_p, R_0, n_exp
		real(kind=AbqRK)::Stress_Vmises
		
		!Variable initializations
		Stress_Vmises = 0.0_AbqRK
		yield_f 	  = 0.0_AbqRK
		R_0 	= mat_properties%ViscoPlastic_prop(1)
		K_p 	= mat_properties%ViscoPlastic_prop(2)
		n_exp 	= mat_properties%ViscoPlastic_prop(3)
		
		!Compute the Von mises eq stress
		Stress_Vmises = compute_VonMises_eq_stress(state_var_trial%Stress,NDI,NTENS)
		
		!Compute the yield criertia
		if(state_var_previous%Plastic_mult < 1e-6_AbqRK) then
			yield_f = (Stress_Vmises/(1-state_var_previous%Damage_var))-K_p*1e-6_AbqRK**n_exp - R_0
		else
			yield_f = (Stress_Vmises/(1-state_var_previous%Damage_var))-K_p*state_var_previous%Plastic_mult**n_exp - R_0
		end if
		
	End Function
	
	Function compute_yield_function_positive(state_var_trial,mat_properties,NTENS,NDI)result(yield_f_positive)
		!=======================================================================
		!Function to compute the value of yield function positive (f +ve)
		!
		!	Inputs:-
		!	Stress_n2			-> updated stress computed with the converged viscoelastic strain value
		!	ViscoPlastic_prop	-> Viscoelastic properties of kelvin branches
		!	NTENS				-> Stress components number
		!	Damage_var_n		-> Damage variable @ n
		!	Plastic_mult_n		-> Plastic multiplier @ n
		!	
		!
		!	Outputs:-
		!   yield_f				-> yield function value
		!=======================================================================
		IMPLICIT NONE
		!Dummy variable declarations
		Type(material_prop_matrix),intent(in)::mat_properties
		Type(state_variables_matrix),intent(in)::state_var_trial
		integer(kind=AbqIK),intent(in)::NTENS,NDI
		real(kind=AbqRK)::yield_f_positive
		
		!Local variables 
		integer(kind=AbqIK)::i,j
		real(kind=AbqRK)::K_p, R_0, n_exp
		real(kind=AbqRK)::Stress_Vmises
		real(kind=AbqRK)::IntermediateValue1
		real(kind=AbqRK)::IntermediateValue2
		real(kind=AbqRK)::Plastic_mult
		real(kind=AbqRK),parameter::TOEL  = 1e-10_AbqRK
		real(kind=AbqRK),parameter::TOEL2  = 1e-6_AbqRK
		
		!Variable initializations
		Stress_Vmises 		 = 0.0_AbqRK
		yield_f_positive 	 = 0.0_AbqRK
		IntermediateValue1	 = 0.0_AbqRK
		IntermediateValue2	 = 0.0_AbqRK
		R_0 				 = mat_properties%ViscoPlastic_prop(1)
		K_p 				 = mat_properties%ViscoPlastic_prop(2)
		n_exp 				 = mat_properties%ViscoPlastic_prop(3)
		
		
		!Compute the Von mises eq stress
		Stress_Vmises = compute_VonMises_eq_stress(state_var_trial%Stress,NDI,NTENS)
		
		!Compute the yield criertia
		if(state_var_trial%Plastic_mult < 1e-6_AbqRK) then
			IntermediateValue1 = (Stress_Vmises/(1-state_var_trial%Damage_var))-K_p*1e-6_AbqRK**n_exp - R_0
		else
			IntermediateValue1 = (Stress_Vmises/(1-state_var_trial%Damage_var))-K_p*state_var_trial%Plastic_mult**n_exp - R_0
		end if
		
		yield_f_positive = (IntermediateValue1 + abs(IntermediateValue1))*0.5_AbqRK
		
	End Function
	
	Function compute_VonMises_eq_stress(Stress,NDI,NTENS)result(Stress_eq)
		!=======================================================================
		!Function to compute Von Mises equivalent stress
		!
		!	Inputs:-
		!	Stress	-> Stress value
		!
		!	Outputs:-
		!   Strain_ve_total	-> Sum of viscoelastic strains from N branches
		!=======================================================================
		IMPLICIT NONE
		!Dummy variable declarations
		integer(kind=AbqIK),intent(in)::NTENS,NDI
		real(kind=AbqRK),intent(in)::Stress(:)
		real(kind=AbqRK)::Stress_eq
		
		!Local variables
		real(kind=AbqRK),dimension(NTENS)::Stress_dev
		real(kind=AbqRK),dimension(6),parameter::Identity_2 = [1.0_AbqRK,1.0_AbqRK,1.0_AbqRK,0.0_AbqRK,0.0_AbqRK,0.0_AbqRK]
		
		!Variable initializations
		Stress_dev		= 0.0_AbqRK
		Stress_eq 	= 0.0_AbqRK
		
		!Compute the deviatoric part
		Stress_dev = Stress - (1.0_AbqRK/3.0_AbqRK)*sum(Stress(1:3))*Identity_2
		
		!Compute the Von mises equivqlent stress
		Stress_eq = sqrt(1.5_AbqRK*dot_product(Stress_dev(1:NDI),Stress_dev(1:NDI))&
							+3.0_AbqRK*dot_product(Stress_dev(NDI+1:NTENS),Stress_dev(NDI+1:NTENS)))
		
	End Function
	
	Function compute_total_energy_density_Y(Strain_total_n2,state_var_trial,mat_properties,NTENS)result(Y_total)
		!=======================================================================
		!Function to Total energy density Y_total = Y_e + Y_vi
		!
		!	Inputs:-
		!	Strain_total_n2		-> Total strain @ n+1
		!	Strain_Ve_n2		-> Visco elastic strain for all K-branches @ n+1
		!	Strain_Vp_n2		-> Viscoplastic strain @ n+1
		!	C_e					-> Elastic stiffness
		!	ViscoElastic_prop	-> Viscoelastic properties
		!	Nu					-> Poisson's ratio
		!
		!	Outputs:-
		!   Y_total				-> Total energy density
		!=======================================================================
		IMPLICIT NONE
		!Dummy variable declarations
		Type(material_prop_matrix),intent(in)::mat_properties
		Type(state_variables_matrix),intent(in)::state_var_trial
		real(kind=AbqRK),intent(in)::Strain_total_n2(:)
		integer(kind=AbqIK),intent(in)::NTENS
		real(kind=AbqRK)::Y_total
		
		
		!Local variables
		real(kind=AbqRK)::Y_e,Y_v 						!Energy density frome lastic and vsicoelastic parts
		real(kind=AbqRK)::Strain_ve_total(NTENS) 		!Total strain from all visco elastic branches
		real(kind=AbqRK)::Strain(NTENS) 			   	!Total strain ( E - sum(E_visco_eleastic) - E_visco_plastic)
		real(kind=AbqRK)::MatrixMulTempHolder(1,NTENS) 	!Temp holder for intermediate result
		real(kind=AbqRK)::Strain_ve_i(NTENS) 			!Viscoelastic strain of i'th kelvin branch
		real(kind=AbqRK)::Cv_i(NTENS,NTENS) 			!Stiffiness of linear spring in i'th Kelvin branch
		integer(kind=AbqIK):: i
		
		!Variable initializations
		Y_e					= 0.0_AbqRK
		Y_v					= 0.0_AbqRK
		Y_total				= 0.0_AbqRK
		Strain_ve_total		= 0.0_AbqRK
		Strain				= 0.0_AbqRK
		MatrixMulTempHolder	= 0.0_AbqRK
		Strain_ve_i			= 0.0_AbqRK
		Cv_i				= 0.0_AbqRK
		
		!Compute total of viscoelastic strain from all branches
		Strain_ve_total = compute_StrainVe_Total(state_var_trial%Strain_Ve,mat_properties%N_k,NTENS)
		
		!Compute tootal strain -> E - Sum(Evi) -Ep
		Strain = Strain_total_n2 - Strain_ve_total -state_var_trial%Strain_Vp
		
		!Compute Ye
		MatrixMulTempHolder = matmul(reshape(Strain,(/1,NTENS/)),mat_properties%C_e)
		Y_e = 0.5_AbqRK*dot_product(reshape(MatrixMulTempHolder,(/NTENS/)),Strain)
		
		!Compute Y_v
		do i=1,mat_properties%N_k
			MatrixMulTempHolder = 0.0_AbqRK
			Strain_ve_i(:) = state_var_trial%Strain_Ve(i,:)
			call assemble_visco_elastic_stiffness(Cv_i,[mat_properties%ViscoElastic_prop(i,1),mat_properties%Nu],NTENS)
			MatrixMulTempHolder = matmul(reshape(Strain_ve_i,(/1,NTENS/)),Cv_i)
			Y_v = Y_v + 0.5_AbqRK*dot_product(reshape(MatrixMulTempHolder,(/NTENS/)),Strain_ve_i)
			
		end do
		
		!Compute total energy density
		Y_total = (Y_e + Y_v )
		
	End Function
	
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
		forall(i=1:6,j=1:6,i==j)Identity_4(i,j)= 1.0_AbqRK
	
	End Function
	
	Function assemble_derivatives_K_matrix(Strain_total_n2,state_var_trial,state_var_previous,mat_properties,DTIME,NDI,N,NTENS)result(K_dd)
		!=======================================================================
		!Subroutine to assemble the deriavtives matrix K_dd
		!
		!	Inputs:-
		!	Strain_total_n2 	-> Total strain @ (n+1)(k+1)
		!	DTIME				-> delta time
		!	N 					-> Number of kelvin branches 	 
		!    
		!	Outputs:-
		!   K_dd				-> Deriavtives @ (n+1)(k+1)
		!=======================================================================
		IMPLICIT NONE
		!Dummy variable declarations
		Type(state_variables_matrix),intent(in)::state_var_previous
		Type(material_prop_matrix),intent(in)::mat_properties
		Type(state_variables_matrix),intent(in)::state_var_trial
		real(kind=AbqRK),intent(in)::Strain_total_n2(:)
		real(kind=AbqRK),intent(in)::DTIME
		integer(kind=AbqIK),intent(in)::NDI,NTENS,N
		real(kind=AbqRK)::K_dd(N*NTENS+NTENS+2,N*NTENS+NTENS+2)
		
		
		!Local variables
		real(kind=AbqRK)::K_ve_i_j(NTENS,NTENS)
		real(kind=AbqRK)::K_vi_p(NTENS,NTENS)
		real(kind=AbqRK)::K_p_vi(NTENS,NTENS)
		real(kind=AbqRK)::K_r_vi(NTENS)
		real(kind=AbqRK)::K_D_vi(NTENS)
		real(kind=AbqRK)::K_p_p(NTENS,NTENS)
		real(kind=AbqRK)::K_p_r(NTENS)
		real(kind=AbqRK)::K_p_D(NTENS)
		real(kind=AbqRK)::K_r_p(NTENS)
		real(kind=AbqRK)::K_D_p(NTENS)
		real(kind=AbqRK)::K_r_r,K_D_r,K_D_D
		real(kind=AbqRK)::Identity_tensor_4(NTENS,NTENS)
		real(kind=AbqRK)::Y_total,S_p,Beta
		integer(kind=AbqIK):: i,j, Start_index,start_index2
		
		!Variable initializations
		K_dd				= 0.0_AbqRK
		K_ve_i_j			= 0.0_AbqRK
		K_vi_p				= 0.0_AbqRK
		K_p_vi				= 0.0_AbqRK
		K_r_vi				= 0.0_AbqRK
		K_D_vi				= 0.0_AbqRK
		K_p_p				= 0.0_AbqRK
		K_p_r 				= 0.0_AbqRK
		K_p_D 				= 0.0_AbqRK
		K_r_p 				= 0.0_AbqRK
		K_D_p 				= 0.0_AbqRK
		Identity_tensor_4	= 0.0_AbqRK
		S_p 				= mat_properties%ViscoPlastic_prop(6)
		Beta 				= mat_properties%ViscoPlastic_prop(7)
		
		!Assemble the sub matrix K_ve to K_dd
		start_index 	= 1
		start_index2 	= 1
		do i=1, N
			do j=1,N
				K_ve_i_j = compute_derivative_K_vi_vj(mat_properties,i,j,DTIME,NTENS)
				K_dd(start_index:,start_index2:) = K_ve_i_j(:,:)
				start_index2 = start_index2+size(K_ve_i_j,2)
				
			end do
				start_index2 = 1
				start_index  = start_index+size(K_ve_i_j,1)
		end do
		
		!Assemble the sub matrices K_vi_p to K_dd
		start_index 	= 1
		start_index2 	= N*NTENS+1
		do i=1, N
		
			K_vi_p = compute_derivative_K_vi_p(mat_properties,i,NTENS)
			K_dd(start_index:,start_index2:) = K_vi_p(:,:)
			
			start_index = start_index+size(K_vi_p,1)
		end do
		
		!Assemble the sub matrices K_vi_r to K_dd
		start_index 	= 1
		start_index2 	= N*NTENS + NTENS + 1
		do i=1, N
		
			K_dd(start_index:,start_index2) = 0.0_AbqRK
			start_index = start_index+NTENS
		end do
		
		
		!Assemble the sub matrices K_vi_D to K_dd
		start_index 	= 1
		start_index2 	= N*NTENS + NTENS + 2
		do i=1, N
		
			K_dd(start_index:,start_index2) = 0.0_AbqRK
			start_index = start_index+NTENS
		end do
		
		!Assemble the sub matrices K_p_vi to K_dd
		start_index 	= N*NTENS+1
		start_index2 	= 1
		K_p_vi = compute_derivative_K_p_vi(state_var_trial,state_var_previous,mat_properties,NDI,NTENS) !Same for all i, so outside loop
		do i=1, N
			K_dd(start_index:,start_index2:) = K_p_vi(:,:)
			start_index2 = start_index2+size(K_p_vi,2)
		end do
		
		!Assemble the sub matrices K_r_vi to K_dd
		start_index 	= N*NTENS+NTENS+1
		start_index2 	= 1
		K_r_vi = compute_derivative_K_r_vi(state_var_trial,mat_properties,DTIME,NDI,NTENS)	 !Same for all i, so outside loop
		do i=1, N
			K_dd(start_index,start_index2:) = K_r_vi(:)
			start_index2 = start_index2+NTENS
		end do
		
		!Assemble the sub matrices K_D_vi to K_dd
		start_index 	= N*NTENS+NTENS+2
		start_index2 	= 1
		do i=1, N
			K_D_vi = compute_derivative_K_D_vi(Strain_total_n2,state_var_trial,state_var_previous,mat_properties,i,N,NTENS)
										
			K_dd(start_index,start_index2:) = K_D_vi(:)
			start_index2 = start_index2+NTENS
		end do
		
		!Assemble the sub matrices K_p_p to K_dd
		start_index 	= N*NTENS+1
		start_index2 	= N*NTENS+1
		Identity_tensor_4 = I_4()
		K_p_p = Identity_tensor_4 + K_p_vi
		K_dd(start_index:,start_index2:) = K_p_p(:,:)
		
		!Assemble the sub matrices K_p_r to K_dd
		start_index 	= N*NTENS+1
		start_index2 	= N*NTENS+NTENS+1
		K_p_r = compute_derivative_K_p_r(state_var_trial,NDI,NTENS)
		K_dd(start_index:,start_index2) = K_p_r(:)
		
		!Assemble the sub matrices K_p_d to K_dd
		start_index 	= N*NTENS+1
		start_index2 	= N*NTENS+NTENS+2
		K_p_D = compute_derivative_K_p_D(state_var_trial,state_var_previous,mat_properties,DTIME,NDI,NTENS)
		K_dd(start_index:,start_index2) = K_p_D(:)
		
		!Assemble the sub matrices K_r_p to K_dd
		start_index 	= N*NTENS+NTENS+1
		start_index2 	= N*NTENS+1
		K_r_p = K_r_vi
		K_dd(start_index,start_index2:) = K_r_p(:)
		
		!Assemble the sub matrices K_r_r to K_dd
		start_index 	= N*NTENS+NTENS+1
		start_index2 	= N*NTENS+NTENS+1
		K_r_r = compute_derivative_K_r_r(state_var_trial,mat_properties,DTIME,NDI,NTENS)
		K_dd(start_index,start_index2) = K_r_r
		
		!Assemble the sub matrices K_r_D to K_dd
		start_index 	= N*NTENS+NTENS+1
		start_index2 	= N*NTENS+NTENS+2
		!K_dd(start_index,start_index2) = 0.0_AbqRK
		K_dd(start_index,start_index2) = compute_derivative_K_r_D(state_var_trial,state_var_previous,mat_properties,DTIME,NDI,NTENS)
		
		!Assemble the sub matrices K_D_p to K_dd
		start_index 	= N*NTENS+NTENS+2
		start_index2 	= N*NTENS+1
		K_D_p = compute_derivative_K_D_p(Strain_total_n2,state_var_trial,state_var_previous,mat_properties,N,NTENS)
		K_dd(start_index,start_index2:) = K_D_p(:)
		
		!Assemble the sub matrices K_D_r to K_dd
		start_index 	= N*NTENS+NTENS+2
		start_index2 	= N*NTENS+NTENS+1
		Y_total = compute_total_energy_density_Y(Strain_total_n2,state_var_trial,mat_properties,NTENS)
		K_D_r = (Y_total/S_p)**Beta
		K_dd(start_index,start_index2) = K_D_r
		
		!Assemble the sub matrices K_D_D to K_dd
		start_index 	= N*NTENS+NTENS+2
		start_index2 	= N*NTENS+NTENS+2
		K_D_D = 2.0_AbqRK * state_var_trial%Damage_var -1.0_AbqRK - state_var_previous%Damage_var
		K_dd(start_index,start_index2) = K_D_D
		
		
	End Function
	
	
	
	Function compute_and_assemble_K_di(Strain_total_n2,state_var_current,state_var_previous,mat_properties,N,NDI,DTIME,NTENS) result(K_di)
		!=======================================================================
		!Function to compute and asseble the derivative matric K_di ( d Phi/ d strain)
		!
		!	Inputs:-
		!	Stress_n2	-> Stress @ n+1
		!	Strain_total_n2	-> Total strain @ n+1
		!	ViscoPlastic_prop	-> Viscoelastic properties
		!
		!	Outputs:-
		!   K_di	-> derivative matrix K_di ( d Phi/ d strain)
		!=======================================================================
		IMPLICIT NONE
		!Dummy variable declarations
		Type(state_variables_matrix),intent(in)::state_var_current
		Type(state_variables_matrix),intent(in)::state_var_previous
		Type(material_prop_matrix),intent(in)::mat_properties
		real(kind=AbqRK),intent(in)::Strain_total_n2(:)
		real(kind=AbqRK),intent(in)::DTIME
		integer(kind=AbqIK),intent(in)::NTENS,N,NDI
		real(kind=AbqRK)::K_di(N*NTENS+NTENS+2,NTENS)
		
		!Local variables
		real(kind=AbqRK)::K_vi_p(NTENS,NTENS)
		real(kind=AbqRK)::K_p_vi(NTENS,NTENS)
		real(kind=AbqRK)::K_r_vi(NTENS)
		real(kind=AbqRK)::K_D_Et(NTENS)
		integer(kind=AbqIK):: i,j,start_index,start_index2
		
		!Variable initializations
		K_di	= 0.0_AbqRK
		K_vi_p = 0.0_AbqRK
		K_p_vi = 0.0_AbqRK
		K_r_vi = 0.0_AbqRK
		K_D_Et = 0.0_AbqRK
		
		!Compute and assemble K_vi_p to K_di
		start_index 	= 1
		do i=1, N
		
			K_vi_p = compute_derivative_K_vi_p(mat_properties,i,NTENS)
			K_di(start_index:,:) = -1.0_AbqRK * K_vi_p(:,:)
			start_index = start_index+size(K_vi_p,1)
		end do
		
		!Compute and assemble K_p_vi to K_di
		start_index 	= N*NTENS+1
		K_p_vi = compute_derivative_K_p_vi(state_var_current,state_var_previous,mat_properties,NDI,NTENS) !Same for all i
		K_di(start_index:,:) = -1.0_AbqRK * K_p_vi(:,:)
		
		!Compute and assemble K_r_vi to K_di
		start_index 	= N*NTENS+NTENS+1
		K_r_vi = compute_derivative_K_r_vi(state_var_current,mat_properties,DTIME,NDI,NTENS)	 !Same for all i
		K_di(start_index,:) = -1.0_AbqRK * K_r_vi(:)
		
		!Compute and assemble K_D_Et
		start_index 	= N*NTENS+NTENS+2
		K_D_Et = -1.0_AbqRK*compute_derivative_K_D_p(Strain_total_n2,state_var_current,state_var_previous,mat_properties,N,NTENS)
		K_di(start_index,:) = K_D_Et(:)
		
		
	End Function
	
	Function compute_and_assemble_K_id(state_var_current,mat_properties,N,NTENS)result(K_id)
		!=======================================================================
		!Function to compute and asseble the derivative matric K_id ( d Stress/ d state variables)
		!
		!	Inputs:-
		!	Stress_n2	-> Stress @ n+1
		!	Plastic_mult_n2 	-> Plastic multiplier @ (n+1)(k+1)
		!	C_e	-> Elastic stiffness
		!
		!	Outputs:-
		!   K_id	-> derivative matrix K_id ( d Stress/ d stte variables)
		!=======================================================================
		IMPLICIT NONE
		!Dummy variable declarations
		Type(state_variables_matrix),intent(in)::state_var_current
		Type(material_prop_matrix),intent(in)::mat_properties
		integer(kind=AbqIK),intent(in)::NTENS,N
		real(kind=AbqRK)::K_id(NTENS,N*NTENS+NTENS+2)
		
		
		!Local variables
		real(kind=AbqRK)::K_sigma_vi(NTENS,NTENS) !Derivative of sigma with viscoealstic strain
		real(kind=AbqRK)::K_sigma_vp(NTENS,NTENS) !Derivative of sigma with viscopalstic strain
		real(kind=AbqRK)::K_sigma_r(NTENS) !Derivative of sigma with plastic multiplier r
		real(kind=AbqRK)::K_sigma_D(NTENS) !Derivative of sigma with plastic multiplier D
		integer(kind=AbqIK):: i,j,start_index,start_index2
		
		!Variable initializations
		K_id	= 0.0_AbqRK
		K_sigma_vi = 0.0_AbqRK
		K_sigma_vp = 0.0_AbqRK
		K_sigma_r = 0.0_AbqRK
		K_sigma_D = 0.0_AbqRK
		
		!Compute and assemble K_sigma_vi to K_di
		start_index 	= 1
		start_index2 	= 1
		do i=1, N
		
			K_sigma_vi = -1.0_AbqRK *(1-state_var_current%Damage_var)*mat_properties%C_e
			K_id(start_index:,start_index2:) = K_sigma_vi(:,:)
			start_index2 = start_index2+size(K_sigma_vi,2)
		end do
	
		
		
		!Compute and assemble K_sigma_vp to K_di
		start_index 	= 1
		start_index2 	= N*NTENS+1
		K_sigma_vp = -1.0_AbqRK *(1-state_var_current%Damage_var)*mat_properties%C_e
		K_id(start_index:,start_index2:) = K_sigma_vp(:,:)
		
		!Compute and assemble K_sigma_r to K_di
		start_index 	= 1
		start_index2 	= N*NTENS+NTENS+1
		K_sigma_r = 0.0_AbqRK
		K_id(start_index:,start_index2) = K_sigma_r(:)
		
		!Compute and assemble K_sigma_D to K_di
		start_index 	= 1
		start_index2 	= N*NTENS+NTENS+2
		K_sigma_D = -1.0_AbqRK *state_var_current%Stress/(1-state_var_current%Damage_var)
		K_id(start_index:,start_index2) = K_sigma_D(:)
		
	End Function
	
	Function compute_projected_stress(state_var_trial,state_var_previous,mat_properties,NDI,NTENS)result(Stress_proj)
		!=======================================================================
		!Function to compute projected stress for stabilising the local NR scheme.
		!
		!	Inputs:-
		!	state_var_trial 			-> trail state
		!	state_var_previous			-> previous converged state
		!
		!	Outputs:-
		!   Stress_proj					-> Projected stress
		!=======================================================================
		IMPLICIT NONE
		!Dummy variable declarations
		Type(state_variables_matrix),intent(in)::state_var_trial,state_var_previous
		Type(material_prop_matrix),intent(in)::mat_properties
		integer(kind=AbqIK),intent(in)::NTENS,NDI
		real(kind=AbqRK)::Stress_proj(NTENS)
		
		!Local variables
		real(kind=AbqRK),dimension(6),parameter::Identity_tensor_2 = [1.0_AbqRK,1.0_AbqRK,1.0_AbqRK,0.0_AbqRK,0.0_AbqRK,0.0_AbqRK]
		real(kind=AbqRK)::Stress_eq, Stress_dev(NTENS), Stress_hyd(NTENS)	!Von Mises equivalent stress and Deviatoric part of Stress
		real(kind=AbqRK)::Stress_proj_dev(NTENS)							!Von Mises equivalent stress and Deviatoric part of Stress
		real(kind=AbqRK)::IntermediateValue1
		real(kind=AbqRK)::H_p,m_exp ,n_exp, K_p, R_0						!Viscoplastic properties (Ref documentation for details)
		
		!Variable initializations
		Stress_proj 			= 0.0_AbqRK
		Stress_proj_dev 		= 0.0_AbqRK
		IntermediateValue1 		= 0.0_AbqRK
		R_0						= mat_properties%ViscoPlastic_prop(1)
		K_p 					= mat_properties%ViscoPlastic_prop(2)
		n_exp 					= mat_properties%ViscoPlastic_prop(3)
		
		!Compute the deviatoric part
		Stress_dev = state_var_trial%Stress - (1.0_AbqRK/3.0_AbqRK)*sum(state_var_trial%Stress(1:3))*Identity_tensor_2
		Stress_hyd = (1.0_AbqRK/3.0_AbqRK)*sum(state_var_trial%Stress(1:3))*Identity_tensor_2
		
		!Compute equivalent stress
		Stress_eq = compute_VonMises_eq_stress(state_var_trial%Stress,NDI,NTENS)
		IntermediateValue1 = Stress_eq/((1-state_var_previous%Damage_var)*(K_p*((state_var_previous%Plastic_mult)**n_exp) + R_0))
		
		!Dev part of projected stress
		Stress_proj_dev = Stress_dev/(IntermediateValue1)
		
		!Stress proj
		Stress_proj = Stress_proj_dev + Stress_hyd
		
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
		if(info .ne. 0) write(6,*) ' Matrix is singular. Not invertable  : '
		
		!If successfull, invert matrix using DGETRI. Uses LU decomposition and ipiv caluclated by DGETRF. Returns LU as the inverted matrix.
		call DGETRI(M,LU,LDA,ipiv,work,LWORK,info)
		if(info .ne. 0) write(6,*) ' Matrix is not inverted  : '
		
		!inverse
		Matrix_inv = LU
		
	End Subroutine
	
	Subroutine linear_system_solver(A,X,B,A_LU)
		!=======================================================================
		!Subroutine to solve a system of linear equations of the form A*X = B.
		!Uses the Intel mkl LAPACK solver : DGESV
		!Ref : Intel® Math Kernel Library Reference Manual
		!
		!	Inputs:-
		!	A	-> Coeffcient matrix
		!	B 	-> RHS
		!
		!	Outputs:-
		!   X 	-> Unknowns
		!=======================================================================
		IMPLICIT NONE
		!Dummy variable declarations
		real(kind=AbqRK),intent(in)::A(:,:)
		real(kind=AbqRK),intent(in)::B(:)
		real(kind=AbqRK),intent(out)::X(size(B,1))
		real(kind=AbqRK),intent(out)::A_LU(size(A,1),size(A,2))
		
		!Local variables
		real(kind=AbqRK)::A_solver(size(A,1),size(A,2))
		real(kind=AbqRK)::B_solver(size(B,1))
		integer(kind=AbqIK)::N ! Order of matrix A
		integer(kind=AbqIK)::NRHS ! Number of columns of matrix B
		integer(kind=AbqIK)::LDA ! Leading dimension of A ie max(1,N)
		integer(kind=AbqIK)::LDB ! Leading dimension of B ie max(1,N)
		integer(kind=AbqIK)::info ! info returned by LAPACK routines
		integer(kind=AbqIK)::ipiv(size(A,1)) !pivot indices from DGETRF. Dimension - min(M,N)
		
		
		!Initializations
		A_solver = A
		B_solver = B
		A_LU = 0.0_AbqRK
		N = size(A,1)
		NRHS =1
		LDA = size(A,1)
		LDB = size(B,1)
		X = 0.0_AbqRK
		
		!Call LAPACK subroutine DGESV to solve for X.
		call DGESV( N, NRHS, A_solver, LDA, IPIV, B_solver, LDB, INFO ) 
		
		if (info .eq. 0) then
			X 	 = B_solver
			A_LU = A_solver
		else
			write(6,*) ' Linear system solver failed to solve for equation A*X = B. From matrix UMAT'
			X = 0.0_AbqRK
		end if 
		
		
	End Subroutine
	
	Subroutine LU_factorization_matrix(A,A_LU,ipiv)
		!=======================================================================
		!Subroutine to solve a system of linear equations of the form A*X = B.
		!Uses the Intel mkl LAPACK solver : DGESV
		!Ref : Intel® Math Kernel Library Reference Manual
		!
		!	Inputs:-
		!	A	-> Coeffcient matrix
		!	B 	-> RHS
		!
		!	Outputs:-
		!   X 	-> Unknowns
		!=======================================================================
		IMPLICIT NONE
		!Dummy variable declarations
		real(kind=AbqRK),intent(inout)::A(:,:)
		real(kind=AbqRK),intent(out)::A_LU(size(A,1),size(A,2))
		integer(kind=AbqIK),intent(out)::ipiv(size(A,1))
		
		!Local variables
		integer(kind=AbqIK)::N ! Order of matrix A
		integer(kind=AbqIK)::M ! Order of matrix A
		integer(kind=AbqIK)::NRHS ! Number of columns of matrix B
		integer(kind=AbqIK)::LDA ! Leading dimension of A ie max(1,N)
		integer(kind=AbqIK)::info ! info returned by LAPACK routines
		
		
		
		!Initializations
		A_LU = 0.0_AbqRK
		M = size(A,1)
		N = size(A,2)
		LDA = size(A,1)
		
		
		!Call LAPACK subroutine DGESV to solve for X.
		call dgetrf( M, N, A, LDA, ipiv, info )
		
		if (info .eq. 0) then
			A_LU = A
		else
			write(6,*) ' Linear system solver failed to LU factorize matrix A. From matrix UMAT'
			A_LU = 0.0_AbqRK
			ipiv = 0.0_AbqRK
		end if 
		
		
	End Subroutine
	


END MODULE Matrix_utilities
