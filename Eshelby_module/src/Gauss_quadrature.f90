!include 'ABQinterface.f90'
!DEC$ FREEFORM
!===============================================================================
! Routine to generate gauss points and weights for a given interval.
! Number of gauss points required and the interval is given as input.
! Routine will return a array of gauss points and an array corresponding weights.
! Ref : Code adapted from Numerical recipies-Art of scientific computing, William H. Press
!===============================================================================

SUBROUTINE gauss_quadrature(interval_x1, interval_x2, gauss_points,weights, num_gauss_points)
	!USE ABQINTERFACE
	IMPLICIT NONE
	!===========================================================================
    !Dummy variable declarations
    !===========================================================================
    INTEGER,PARAMETER::AbqRK=8,AbqIK=4
    real(kind=AbqRK)::interval_x1,interval_x2
    integer(kind=AbqIK)::num_gauss_points
    real(kind=AbqRK),dimension(0:num_gauss_points-1)::gauss_points,weights
    !f2py intent(in) interval_x1, interval_x2,num_gauss_points
	!f2py intent(out) gauss_points,weights
	!f2py depend(num_gauss_points)gauss_points
    !f2py depend(num_gauss_points)weights
   
    
    !===========================================================================
    !Local variable declarations
    !===========================================================================
     integer(kind=AbqIK) :: half_size,loop_index_i,loop_index_j
     real(kind=AbqRK)::lower_bound,upper_bound,z,z1,tolerence,temp,zeta_1_func
     real(kind=AbqRK)::pp,p3,p2,p1 !p1 desired legendre polynomial, pp its derivative
  
    !Variable initializations
    gauss_points	= 0.0_AbqRK
    weights			= 0.0_AbqRK
    half_size		= (num_gauss_points+1)/2
    upper_bound		= 0.5_AbqRK *(interval_x2+interval_x1)
    lower_bound		= 0.5_AbqRK *(interval_x2-interval_x1)
    tolerence		= 1.0e-14
    z = 0.0_AbqRK
    z1 = 0.0_AbqRK
    p3 = 0.0_AbqRK  
    !Start of routine logic
    DO loop_index_i = 0, half_size-1
		z = COS(3.141592654_AbqRK*(loop_index_i+0.75_AbqRK)/(num_gauss_points+0.5_AbqRK))
		
		DO
			p1 = 1.0_AbqRK
			p2 = 0.0_AbqRK
			DO loop_index_j = 0, num_gauss_points-1
				p3 = p2
				p2 = p1
				p1 = ((2.0_AbqRK*loop_index_j+1.0_AbqRK)*z*p2-loop_index_j*p3)/(loop_index_j+1)
			END DO
			pp = num_gauss_points*(z*p1-p2)/(z*z-1.0_AbqRK)
			z1 = z
			z = z1-(p1/pp)
			IF (ABS(z-z1) <tolerence) then
				EXIT
			END IF
		END DO
		gauss_points(int(loop_index_i))=upper_bound-lower_bound*z
		gauss_points(int(num_gauss_points-1-loop_index_i)) = upper_bound+lower_bound*z
		weights(int(loop_index_i)) = 2.0_AbqRK*lower_bound/((1.0_AbqRK - z*z)*pp*pp)
		weights(int(num_gauss_points-1-loop_index_i)) = weights(int(loop_index_i))
		
    END DO
	
END SUBROUTINE gauss_quadrature
