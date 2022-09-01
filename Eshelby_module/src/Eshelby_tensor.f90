!DEC$ FREEFORM
include 'Gauss_quadrature.f90'
include 'Eshelby_utilities.f90'
SUBROUTINE Eshelby_tensor(C_0,a,S_E,M_gauss,N_gauss)
	!===========================================================================
	! Routine to calculate components of Eshelby tensor. Returns SE in Voigt form.
	! Inputs :
	! C_0 - Stiffness of virgin material.
	! a - Principal dimensions of the Inclusion (a1,a2,a3)
	! M_gauss - Number of gauss points for variable zeta3
	! N_gauss - Number of gauss points for variable omega
	!===========================================================================
	!USE ABQINTERFACE
	use utilities
	IMPLICIT NONE
	!===========================================================================
    !Dummy variable declarations
    !===========================================================================
    INTEGER,PARAMETER::AbqRK=8,AbqIK=4
    !INTEGER,PARAMETER::AbqRK=KIND(8),AbqIK=KIND(4)
    real(Kind=AbqRK),intent(in) ::C_0(6,6)
    real(Kind=AbqRK),intent(in)::a(3)
    integer(Kind=AbqIK),intent(in):: M_gauss,N_gauss
    real(Kind=AbqRK),intent(out)::S_E(6,6)
    !f2py intent(in) C_0,a,M_gauss,N_gauss
	!f2py intent(out) S_E
    !===========================================================================
    !Local variable declarations
    !===========================================================================
    integer(Kind=AbqIK):: index_1,index_2
    integer(Kind=AbqIK):: i,j,k,l
    real(Kind=AbqRK)::weights_zeta(M_gauss),weights_omega(N_gauss)
    real(Kind=AbqRK)::gauss_points_zeta(M_gauss),gauss_points_omega(N_gauss)
    real(kind=AbqRK)::interval_x1_zeta,interval_x2_zeta
    real(kind=AbqRK)::interval_x1_omega,interval_x2_omega
    real(Kind=AbqRK),parameter::pi=4*ATAN(1.d0)
    !Variable initializations
    S_E = 0.0d0
    gauss_points_zeta = 0.0d0
    gauss_points_omega = 0.0d0
    weights_zeta = 0.0d0
    weights_omega = 0.0d0
    !pi = 3.141592654_AbqRK
    interval_x1_omega = 0.0d0
    interval_x2_omega = 2*pi
    interval_x1_zeta = -1.0d0
    interval_x2_zeta = 1.0d0
    !Get weights and gauss points for zeta3
    call gauss_quadrature(interval_x1_zeta, interval_x2_zeta,gauss_points_zeta,weights_zeta,M_gauss)
    !Get weights and gauss points for omega
    call gauss_quadrature(interval_x1_omega, interval_x2_omega,gauss_points_omega,weights_omega,N_gauss)
    
    DO index_1 = 1,6
		DO index_2 = 1,6
			i = 0
			j = 0
			k = 0
			l = 0
			! Get the tensorial index corresponding to the Voigt index used. All calculations to find SE_ijkl are done using tensorial operations.
			call Return_Tensorial_index(i,j,k,l,index_1,index_2)
			!Call the function that calculates the Eshelby tensor componenst SE_i,j,k,l
			S_E(index_1,index_2)= SE_ijkl_func(i,j,k,l,C_0,a,weights_zeta,&
					weights_omega,gauss_points_zeta,gauss_points_omega,M_gauss,N_gauss)
		END DO
    END DO
    !Mutiplying the weight of 2 to the lower half of Eshelby tensor rep in voigt notation
    S_E(4,:)= 2*S_E(4,:)
    S_E(5,:)= 2*S_E(5,:)
    S_E(6,:)= 2*S_E(6,:)
   
    
END SUBROUTINE Eshelby_tensor
