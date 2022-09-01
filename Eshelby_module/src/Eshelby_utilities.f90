!DEC$ FREEFORM
MODULE utilities
	implicit none
	contains
	!---------------------------------------------------------------------------
	FUNCTION zeta_1_func(omega,zeta_3)result(zeta_1)
		!===========================================================================
		! Function to calculate Zeta_1.
		!===========================================================================
		!USE ABQINTERFACE
		IMPLICIT NONE
		!===========================================================================
		!Dummy variable declarations
		!===========================================================================
		INTEGER,PARAMETER::AbqRK=8,AbqIK=4
		real(Kind=AbqRK),intent(in) ::omega,zeta_3
		real(Kind=AbqRK)::zeta_1
		!f2py intent(in) omega,zeta_3
		!f2py intent(out) zeta_1
		!Variable initializations
		zeta_1 = 0.0d0
		
		zeta_1 = SQRT(1- zeta_3**2)*COS(omega)
		
	END FUNCTION
	
	FUNCTION zeta_2_func(omega,zeta_3)result(zeta_2)
		!===========================================================================
		! Function to calculate Zeta_2.
		!===========================================================================
		!USE ABQINTERFACE
		IMPLICIT NONE
		!===========================================================================
		!Dummy variable declarations
		!===========================================================================
		INTEGER,PARAMETER::AbqRK=8,AbqIK=4
		real(Kind=AbqRK) ::omega
		real(Kind=AbqRK)::zeta_3,zeta_2
		!f2py intent(in) omega,zeta_3
		!f2py intent(out) zeta_2
		!Variable initializations
		zeta_2 = 0.0d0
		
		zeta_2 = SQRT(1- zeta_3**2)*SIN(omega)
	END FUNCTION zeta_2_func
	
	FUNCTION zeta_bar_func(omega,zeta_3,a,i) result(zeta_bar)
		!===========================================================================
		! Function to calculate Zeta_bar.
		!===========================================================================
		!USE ABQINTERFACE
		IMPLICIT NONE
		!===========================================================================
		!Dummy variable declarations
		!===========================================================================
		INTEGER,PARAMETER::AbqRK=8,AbqIK=4
		real(Kind=AbqRK) ::omega
		real(Kind=AbqRK)::zeta_3,zeta_bar
		real(Kind=AbqRK) ::a(3)
		integer(Kind=AbqIK)::i
		!f2py intent(in) omega,zeta_3,a,i
		!f2py intent(out) zeta_bar
		!Variable initializations
		zeta_bar = 0.0d0
		
		select case(i)
			case(1)
				zeta_bar= zeta_1_func(omega,zeta_3)/a(1)
			case(2)
				zeta_bar= zeta_2_func(omega,zeta_3)/a(2)
			case(3)
				zeta_bar= zeta_3/a(3)
		end select
		
	END FUNCTION zeta_bar_func
	
	FUNCTION Voigt_map_func(index_part)result(Voigt_map)
		!===========================================================================
		! Helper function for the routine Return_Voigt_index. Defines the mapping
		! between tensor indices and matrix indices.
		!===========================================================================
		!USE ABQINTERFACE
		IMPLICIT NONE
		!===========================================================================
		!Dummy variable declarations
		!===========================================================================
		INTEGER,PARAMETER::AbqRK=8,AbqIK=4
		integer(Kind=AbqIK)::Voigt_map,index_part
		!f2py intent(in) index_part
		!f2py intent(out) Voigt_map
		!Variable initializations
		Voigt_map = 0
		select case(index_part)
			case(11)
				Voigt_map= 1
			case(22)
				Voigt_map= 2
			case(33)
				Voigt_map= 3
			case(12,21)
				Voigt_map= 4
			case(13,31)
				Voigt_map= 5
			case(23,32)
				Voigt_map= 6
		end select
	END FUNCTION Voigt_map_func
	
	SUBROUTINE Return_Voigt_index(i,j,k,l,index_1,index_2)
		!===========================================================================
		! This routine returns the voigt index corresponding to an index of a fourth
		! tensor. For eg C_1111 -> C_11, C_2312 -> C_46.
		! This is used to get the stifness values from C_0, which is an input given in
		! Voigt form. 
		! i,j,k,l : are the index of the 4 th order tensor component
		! index_1,index_2 : are the indices of the corresponding matrix in Voigt form
		!===========================================================================
		!USE ABQINTERFACE
		IMPLICIT NONE
		!===========================================================================
		!Dummy variable declarations
		!===========================================================================
		INTEGER,PARAMETER::AbqRK=8,AbqIK=4
		integer(Kind=AbqIK),intent(in)::i,j,k,l
		integer(Kind=AbqIK),intent(out)::index_1,index_2
		!f2py intent(in) i,j,k,l
		!f2py intent(out) index_1,index_2
		!===========================================================================
		!Local variable declarations
		!===========================================================================
		integer(Kind=AbqIK)::part1,part2
		!Variable initializations
		part1= 0
		part2= 0
		index_1=0
		index_1=0
		
		
		part1= i*10+j
		part2= k*10+l
		index_1 = Voigt_map_func(part1)
		index_2 = Voigt_map_func(part2)
		
	END SUBROUTINE Return_Voigt_index
	
	FUNCTION K_ik_func(C_0,i,k,omega,zeta_3,a)result(K_ik)
		!===========================================================================
		! Function to calculate K_ik.
		!===========================================================================

		!USE ABQINTERFACE
		IMPLICIT NONE
		!===========================================================================
		!Dummy variable declarations
		!===========================================================================
		INTEGER,PARAMETER::AbqRK=8,AbqIK=4
		real(Kind=AbqRK) ::C_0(6,6)
		integer(Kind=AbqIK)::i,k
		real(Kind=AbqRK)::K_ik
		real(Kind=AbqRK) ::omega
		real(Kind=AbqRK)::zeta_3
		real(Kind=AbqRK) ::a(3)
		!f2py intent(in) C_0,i,k,omega,zeta_3,a
		!f2py intent(out) K_ik
		!===========================================================================
		!Local variable declarations
		!===========================================================================
		integer(Kind=AbqIK)::j,l
		integer(Kind=AbqIK)::index_1,index_2
		!Variable initializations
		K_ik = 0.0d0
		
		
		DO j=1,3
			DO l=1,3
				index_1 = 0
				index_2 = 0
				call Return_Voigt_index(i,j,k,l,index_1,index_2)
				K_ik = K_ik + C_0(index_1,index_2)*zeta_bar_func(omega,zeta_3,a,j)&
						*zeta_bar_func(omega,zeta_3,a,l)
			END DO
		END DO
		
	END FUNCTION K_ik_func
	
	FUNCTION Permutation_tensor_func(i,j,k) result(Permutation_tensor)
		!===========================================================================
		! Function return the value of permuation tensor for a given set of indicies.
		!===========================================================================
		!USE ABQINTERFACE
		IMPLICIT NONE
		!===========================================================================
		!Dummy variable declarations
		!===========================================================================
		INTEGER,PARAMETER::AbqRK=8,AbqIK=4
		integer(Kind=AbqIK)::i,j,k,Permutation_tensor
		!f2py intent(in) i,j,k
		!f2py intent(out) Permutation_tensor
		!===========================================================================
		!Local variable declarations
		!===========================================================================
		integer(Kind=AbqIK)::case_non_zero
		!Variable initializations
		Permutation_tensor=0
		
		if ( i==j .or. j==k .or. k==i) then
			Permutation_tensor = 0
		else
			case_non_zero = (100*i)+(10*j)+k
			select case(case_non_zero)
				case(123,231,312)
					Permutation_tensor= 1
				case(132,321,213)
					Permutation_tensor= -1
			end select
		end if
		
	END FUNCTION Permutation_tensor_func
	
	FUNCTION N_ij_func(i,j,C_0,omega,zeta_3,a)result(N_ij)
		!===========================================================================
		! Function calculates the value of N_ij
		!===========================================================================
		!USE ABQINTERFACE
		IMPLICIT NONE
		!===========================================================================
		!Dummy variable declarations
		!===========================================================================
		INTEGER,PARAMETER::AbqRK=8,AbqIK=4
		integer(Kind=AbqIK)::i,j
		real(Kind=AbqRK)::C_0(6,6)
		real(Kind=AbqRK)::N_ij
		real(Kind=AbqRK) ::omega
		real(Kind=AbqRK)::zeta_3
		real(Kind=AbqRK) ::a(3)
		!f2py intent(in) i,j,C_0,omega,zeta_3,a
		!f2py intent(out) N_ij
		!===========================================================================
		!Local variable declarations
		!===========================================================================
		integer(Kind=AbqIK)::k,l,m,n
		!Variable initializations
		N_ij=0.0d0
		
		DO k=1,3
			DO l=1,3
				DO m=1,3
					DO n=1,3
					N_ij = N_ij + (0.5_AbqRK*Permutation_tensor_func(i,k,l)*&
							Permutation_tensor_func(j,m,n)*K_ik_func(C_0,k,m,omega,zeta_3,a)*&
							K_ik_func(C_0,l,n,omega,zeta_3,a))
					END DO
				END DO
			END DO
		END DO

	END FUNCTION N_ij_func
	
	FUNCTION D_func(C_0,omega,zeta_3,a) result(D)
		!===========================================================================
		! Function calculates the value of D
		!===========================================================================
		!USE ABQINTERFACE
		IMPLICIT NONE
		!===========================================================================
		!Dummy variable declarations
		!===========================================================================
		INTEGER,PARAMETER::AbqRK=8,AbqIK=4
		real(Kind=AbqRK)::C_0(6,6)
		real(Kind=AbqRK)::D
		real(Kind=AbqRK) ::omega
		real(Kind=AbqRK)::zeta_3
		real(Kind=AbqRK) ::a(3)
		!f2py intent(in) C_0,omega,zeta_3,a
		!f2py intent(out) D
		!===========================================================================
		!Local variable declarations
		!===========================================================================
		integer(Kind=AbqIK)::m,n,l
		!Variable initializations
		D=0.0d0
		
		DO m=1,3
			DO n=1,3
				DO l=1,3
					D = D + (Permutation_tensor_func(m,n,l)*&
							 K_ik_func(C_0,m,1,omega,zeta_3,a)*&
							 K_ik_func(C_0,n,2,omega,zeta_3,a)*&
							 K_ik_func(C_0,l,3,omega,zeta_3,a))
				END DO
			END DO
		END DO
		
	END FUNCTION D_func

	FUNCTION G_func(i,j,k,l,C_0,omega,zeta_3,a) result(G)
		!===========================================================================
		! Function calculates the value of G
		!===========================================================================
		!USE ABQINTERFACE
		IMPLICIT NONE
		!===========================================================================
		!Dummy variable declarations
		!===========================================================================
		INTEGER,PARAMETER::AbqRK=8,AbqIK=4
		integer(Kind=AbqIK)::i,j,k,l
		real(Kind=AbqRK)::C_0(6,6)
		real(Kind=AbqRK)::G
		real(Kind=AbqRK) ::omega
		real(Kind=AbqRK)::zeta_3
		real(Kind=AbqRK) ::a(3)
		!f2py intent(in) i,j,k,l,C_0,omega,zeta_3,a
		!f2py intent(out) G
		!===========================================================================
		!Local variable declarations
		!===========================================================================
		!Variable initializations
		G=0.0d0
		
		G = (zeta_bar_func(omega,zeta_3,a,k)*zeta_bar_func(omega,zeta_3,a,l)*&
				N_ij_func(i,j,C_0,omega,zeta_3,a))/D_func(C_0,omega,zeta_3,a)
		
	END FUNCTION G_func
	
	SUBROUTINE Voigt_map_reverse(index_1,i,j)
		!===========================================================================
		! Helper function for the routine Return_Tensorial_index. Defines the mapping
		! between tensor indices and matrix indices.
		!===========================================================================
		
		!USE ABQINTERFACE
		IMPLICIT NONE
		!===========================================================================
		!Dummy variable declarations
		!===========================================================================
		INTEGER,PARAMETER::AbqRK=8,AbqIK=4
		integer(Kind=AbqIK),intent(in)::index_1
		integer(Kind=AbqIK),intent(out)::i,j
		!f2py intent(in) index_1
		!f2py intent(out) i,j
		!Variable initializations
		
		select case(index_1)
			case(1)
				i=1
				j=1
			case(2)
				i=2
				j=2
			case(3)
				i=3
				j=3
			case(4)
				i=1
				j=2
			case(5)
				i=1
				j=3
			case(6)
				i=2
				j=3
		end select
	END SUBROUTINE Voigt_map_reverse
	
	SUBROUTINE Return_Tensorial_index(i,j,k,l,index_1,index_2)
		!===========================================================================
		! This routine returns the index of a fourth order tensor corresponding to a voigt
		! index  . For eg C_11 -> C_1111, C_46 -> C_2312.
		! This is used in the calculation of SE_ijkl. Eshelby tensor is returned by the
		! subroutine in Voigt form. Only the 36 components are calculated, hence this transformation is req. 
		! i,j,k,l : are the index of the 4 th order tensor component
		! index_1,index_2 : are the indices of the corresponding matrix in Voigt form
		!===========================================================================
		!USE ABQINTERFACE
		IMPLICIT NONE
		!===========================================================================
		!Dummy variable declarations
		!===========================================================================
		INTEGER,PARAMETER::AbqRK=8,AbqIK=4
		integer(Kind=AbqIK),intent(in)::index_1,index_2
		integer(Kind=AbqIK),intent(out)::i,j,k,l
		!f2py intent(in) index_1,index_2
		!f2py intent(out) i,j,k,l
		!===========================================================================
		!Local variable declarations
		!===========================================================================
		integer(Kind=AbqIK)::part1,part2
		!Variable initializations
		i = 0
		j = 0
		k = 0
		l = 0
		call Voigt_map_reverse(index_1,i,j)
		call Voigt_map_reverse(index_2,k,l)
		
	END SUBROUTINE Return_Tensorial_index
	
	FUNCTION SE_ijkl_func(i,j,k,l,C_0,a,weights_zeta,weights_omega,gauss_points_zeta,&
						gauss_points_omega,M_gauss,N_gauss) result(SE_ijkl)
	!USE ABQINTERFACE
	IMPLICIT NONE
	!===========================================================================
	!Dummy variable declarations
	!===========================================================================
	INTEGER,PARAMETER::AbqRK=8,AbqIK=4
	integer(Kind=AbqIK):: i,j,k,l
	integer(Kind=AbqIK):: M_gauss,N_gauss,index_1,index_2
	real(Kind=AbqRK)::weights_zeta(M_gauss),weights_omega(N_gauss)
	real(Kind=AbqRK)::gauss_points_zeta(M_gauss),gauss_points_omega(N_gauss)
	real(Kind=AbqRK) ::C_0(6,6),SE_ijkl
	real(Kind=AbqRK)::a(3)
	!f2py intent(in) i,j,k,l,weights_zeta,weights_omega,gauss_points_zeta
	!f2py intent(in) weights_zeta,weights_omega,gauss_points_zeta
	!f2py intent(in) gauss_points_omega,M_gauss,N_gauss,C_0,a
	!f2py intent(out) SE_ijkl
	!f2py depend(M_gauss)weights_zeta
    !f2py depend(N_gauss)weights_omega
    !f2py depend(M_gauss)gauss_points_zeta
    !f2py depend(N_gauss)gauss_points_omega
	!===========================================================================
	!Local variable declarations
	!===========================================================================
	integer(Kind=AbqIK):: m,n,p,q
	real(Kind=AbqRK),parameter::pi=4*ATAN(1.d0)
	real(Kind=AbqRK)::G1,G2
	!Variable initializations
	SE_ijkl = 0.0d0
	index_1 =0
	index_2 =0
	!pi = 3.141592654_AbqRK
	!Perfrom the  numerical integration
	DO m=1,3
		DO n=1,3
			!Fetch the Voigt index corresponding to the tensorial index m,n,k,l to get the Stifness quantity from C_0 which is given in Voigt form.
			index_1 =0
			index_2 =0
			call Return_Voigt_index(m,n,k,l,index_1,index_2)
			DO p=1,M_gauss
				DO q=1,N_gauss
				G1 = 0.0d0
				G2 = 0.0d0
				G1 = G_func(i,m,j,n,C_0,gauss_points_omega(q),gauss_points_zeta(p),a)
				G2 = G_func(j,m,i,n,C_0,gauss_points_omega(q),gauss_points_zeta(p),a)
				!Compute SE_ijkl Equation ref : Hybrid micromechanical-phenomenological modelling... Apendix A
				SE_ijkl = SE_ijkl + ((1.0d0/(8*pi))*C_0(index_1,index_2)* &
						  ( G1 + G2)*weights_omega(q)*weights_zeta(p))
				END DO
			END DO
		END DO
	END DO
		
	END FUNCTION SE_ijkl_func
	
	
END MODULE
