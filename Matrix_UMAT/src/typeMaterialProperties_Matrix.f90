!DEC$ FREEFORM
!===============================================================================
! Class(type) definition for the matrix material properties. Member functions to 
! allocate and deallocate the varibles are included.
!===============================================================================
Module Type_material_prop_matrix
	USE ABQINTERFACE
    IMPLICIT NONE
	Public
	
	Type material_prop_matrix
		real(kind=AbqRK),dimension(:,:),allocatable::C_e				!Elastic stifness in Viogt notation
		real(kind=AbqRK),dimension(:,:),allocatable::ViscoElastic_prop	!Viscoelastic material properties
		real(kind=AbqRK),dimension(:),allocatable::ViscoPlastic_prop	!Viscoplastic and damage properties
		real(kind=AbqRK)::E_e,Nu 										!Youngs Modulus of linear spring an poisson's ratio
		integer(kind=AbqIK) :: N_k
		
		Contains
		
		Procedure::allocation 	=> allocate_material_prop_matrix
		Procedure::deallocation => deallocate_material_prop_matrix
		Procedure::get_elastic_stiffness => assemble_elastic_stiffness
	
	end Type material_prop_matrix
			  
	
	Private ::  allocate_material_prop_matrix,deallocate_material_prop_matrix
			  
	Contains
	
	Subroutine allocate_material_prop_matrix(this,PROPS,NTENS)
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
		Class(material_prop_matrix),intent(inout)::this
		real(kind=AbqRK),intent(in)::PROPS(:)
		integer(kind=AbqIK),intent(in)::NTENS
		
		!Local variables
		integer(kind=AbqIK) ::i
		
		!Get number of Kelvin branches
		this%N_k = int(PROPS(3))
		
		!Allocate the varibles
		allocate(this%C_e(NTENS,NTENS))
		allocate(this%ViscoElastic_prop(this%N_k,2))
		allocate(this%ViscoPlastic_prop(7))
		
		!Initialize the varibles
		this%C_e 				= 0.0_AbqRK
		this%ViscoElastic_prop 	= 0.0_AbqRK
		this%ViscoPlastic_prop 	= 0.0_AbqRK
		
		!Set Youngs modulus and poisson's ratio
		this%E_e = PROPS(1)
		this%Nu  = PROPS(2)
		
		!Assign viscoelastic and viscoplastic properties
		forall(i=1:this%N_k)this%ViscoElastic_prop(i,:)= PROPS((2*i+2):)
		this%ViscoPlastic_prop = PROPS(3+this%N_k*2+1: 10+this%N_k*2)
		
		!Assemble the elastic stiffness from material props
		call this%get_elastic_stiffness(PROPS,NTENS)
		
	End Subroutine
	
	Subroutine deallocate_material_prop_matrix(this)
		!=======================================================================
		!Subroutine to de-allocate the state variables
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
		Class(material_prop_matrix),intent(inout)::this
	
		
		!Allocate the varibles
		deallocate(this%C_e)
		deallocate(this%ViscoElastic_prop)
		deallocate(this%ViscoPlastic_prop)
		
	End Subroutine
	
	Subroutine assemble_elastic_stiffness(this,PROPS,NTENS)
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
		Class(material_prop_matrix),intent(inout)::this
		real(kind=AbqRK),intent(in)::PROPS(:)
		integer(kind=AbqIK),intent(in)::NTENS
		
		!Local variables
		real(kind=AbqRK)::Lambda,mu
		real(kind=AbqRK)::E_e,nu
		integer(kind=AbqIK)::i,j
		
		!Variable initialization
		this%C_e = 0.0_AbqRK
		
		!Compute Lame's constants : Lambda and mu
		Lambda =  this%E_e*this%Nu/((1-2*this%Nu)*(1+this%Nu))
		mu = this%E_e/(2*(1+this%Nu))
		
		!First half of diagonals
		forall(i=1:3)this%C_e(i,i)= Lambda +2*mu
		!Second half of diagonals
		forall(i=4:6)this%C_e(i,i)= mu
		!Off diagonal terms of upper matrix
		forall(i=1:3,j=1:3,i/=j)this%C_e(i,j)= Lambda
	
	End Subroutine
end Module Type_material_prop_matrix
