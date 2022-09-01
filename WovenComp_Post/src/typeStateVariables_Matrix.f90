!DEC$ FREEFORM
!===============================================================================
! Class(type) definition for the matrix state varibles. Member functions to 
! allocate and deallocate the varibles are included.
!===============================================================================
Module Type_State_variables_matrix
	USE ABQINTERFACE
    IMPLICIT NONE
	Public
	
	Type state_variables_matrix
		real(kind=AbqRK),dimension(:),allocatable::Stress 		!Stress
		real(kind=AbqRK),dimension(:,:),allocatable::Strain_Ve 	!Viscoelastic strain 
		real(kind=AbqRK),dimension(:),allocatable::Strain_Vp	!Viscoplastic strain
		real(kind=AbqRK)::Plastic_mult							!Plastic multiplier 
		real(kind=AbqRK)::Damage_var							!Damage variable
		
		Contains
		
		Procedure::allocation 	=> allocate_state_variables
		Procedure::deallocation => deallocate_state_variables
	
	end Type state_variables_matrix
			  
	
	Private :: allocate_state_variables, deallocate_state_variables  
			  
	Contains
	
	Subroutine allocate_state_variables(this,NTENS,N)
		!=======================================================================
		!Subroutine to allocate the state variables
		!PROPS
		!	Inputs:-
		!	NTENS	-> Number of stress components
		!	N 		-> Number of Kelvin branches
		!              
		!	Outputs:-
		!   Allocates the variables
		!=======================================================================
		IMPLICIT NONE
		!Dummy variable declarations
		Class(state_variables_matrix),intent(inout)::this
		integer(kind=AbqIK),intent(in)::NTENS, N
		
		!Allocate the varibles
		allocate(this%Stress(NTENS))
		allocate(this%Strain_Ve(N,NTENS))
		allocate(this%Strain_Vp(NTENS))
		
		!Initialize the vraibles
		this%Stress 		= 0.0_AbqRK
		this%Strain_Ve 		= 0.0_AbqRK
		this%Strain_Vp 		= 0.0_AbqRK
		this%Plastic_mult 	= 0.0_AbqRK
		this%Damage_var 	= 0.0_AbqRK
		
	End Subroutine
	
	Subroutine deallocate_state_variables(this)
		!=======================================================================
		!Subroutine to deallocate the state variables
		!PROPS
		!	Inputs:-
		!	this	-> Current instance variable
		!              
		!	Outputs:-
		!   Allocates the variables
		!=======================================================================
		IMPLICIT NONE
		!Dummy variable declarations
		Class(state_variables_matrix),intent(inout)::this
		
		!Allocate the varibles
		deallocate(this%Stress)
		deallocate(this%Strain_Ve)
		deallocate(this%Strain_Vp)
		
	End Subroutine

end Module Type_State_variables_matrix
