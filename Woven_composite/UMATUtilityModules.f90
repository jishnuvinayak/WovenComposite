!DEC$ FREEFORM
!===============================================================================
! Seperate modules to include the material routines for the yarns and the matrix of
! woven composite.
!===============================================================================
!Include the lapapck lib used in matrix material routine.
!Include 'lapack.f90'
include 'typeMaterialProperties_Matrix.f90'
include 'typeStateVariables_Matrix.f90'
include 'UMATYarn_Utilities.f90'
include 'UMATMatrix_Utilities.f90'
!===============================================================================
!Module for yarn
!===============================================================================
MODULE Yarn
	USE ABQINTERFACE
	USE Yarn_utilities, only: invert_matrix
    IMPLICIT NONE
    Private
        Public :: UMAT,compute_rotation_matrix,rotate_fourth_order_tensor
        Contains
        !Include yarn material routine.
        Include 'UMATYarn.f90'
     
     Function compute_rotation_matrix(e_yarn_1,e_yarn_2) result(rotation_matrix)
		!=======================================================================
		!Function to compute the rotation matrix using yarn orientation.
		!
		!	Inputs:-
		!	e_yarn_1			-> basis vector along yarn direction
		!	e_yarn_2 			-> basis vector perpendicular to yarn direction
		!	
		!
		!	Outputs:-
		!   rotation_matrix		-> Rotation matrix
		!=======================================================================
		IMPLICIT NONE
		!Dummy variable declarations
		real(kind=AbqRK),intent(in)::e_yarn_1(:)
		real(kind=AbqRK),intent(in)::e_yarn_2(:)
		real(kind=AbqRK)::rotation_matrix(3,3)
		
		!Local variable declarations
		real(kind=AbqRK)::e_yarn_3(3) 			!Third basis vector in yarn local basis system
		real(kind=AbqRK)::e_yarn_1_norm(3),e_yarn_2_norm(3),e_yarn_3_norm(3)
		real(kind=AbqRK)::e_ref_1(3),e_ref_2(3),e_ref_3(3)
		integer(kind=AbqIK)::i,j 				!Running indices
		
		!Variable initializations
		e_yarn_3			= 0.0_AbqRK
		rotation_matrix		= 0.0_AbqRK
		e_ref_1	= [1.0_AbqRK,0.0_AbqRK,0.0_AbqRK]
		e_ref_2	= [0.0_AbqRK,1.0_AbqRK,0.0_AbqRK]
		e_ref_3	= [0.0_AbqRK,0.0_AbqRK,1.0_AbqRK]
		
		!Normalize
		e_yarn_1_norm(:)=e_yarn_1(:)/(sqrt(e_yarn_1(1)**2+e_yarn_1(2)**2+e_yarn_1(3)**2))
		e_yarn_2_norm(:)=e_yarn_2(:)/(sqrt(e_yarn_2(1)**2+e_yarn_2(2)**2+e_yarn_2(3)**2))
		
		!Third basis vector in yarn local cordiante system e3* = e1* x e2*
		e_yarn_3(1) = (e_yarn_1_norm(2)*e_yarn_2_norm(3) - e_yarn_1_norm(3)*e_yarn_2_norm(2) )
		e_yarn_3(2) = (e_yarn_1_norm(3)*e_yarn_2_norm(1) - e_yarn_1_norm(1)*e_yarn_2_norm(3) )
		e_yarn_3(3) = (e_yarn_1_norm(1)*e_yarn_2_norm(2) - e_yarn_1_norm(2)*e_yarn_2_norm(1) )
		
		!Normalize
		e_yarn_3_norm(:)=e_yarn_3(:)/(sqrt(e_yarn_3(1)**2+e_yarn_3(2)**2+e_yarn_3(3)**2))
		
		
		!Compute the roation matrix R_ij = e_ref_i . e_yarn_j
		
		rotation_matrix(1,1) = dot_product(e_ref_1,e_yarn_1_norm)
		rotation_matrix(1,2) = dot_product(e_ref_2,e_yarn_1_norm)
		rotation_matrix(1,3) = dot_product(e_ref_3,e_yarn_1_norm)
		
		rotation_matrix(2,1) = dot_product(e_ref_1,e_yarn_2_norm)
		rotation_matrix(2,2) = dot_product(e_ref_2,e_yarn_2_norm)
		rotation_matrix(2,3) = dot_product(e_ref_3,e_yarn_2_norm)
		
		rotation_matrix(3,1) = dot_product(e_ref_1,e_yarn_3_norm)
		rotation_matrix(3,2) = dot_product(e_ref_2,e_yarn_3_norm)
		rotation_matrix(3,3) = dot_product(e_ref_3,e_yarn_3_norm)
		
	
	End Function
     
     Pure Function mod_r(k,l) result(index_1)
            !=======================================================================
            !Function to compute the rotation matrix using yarn orientation.
            !
            !	Inputs:-
            !	e_yarn_1			-> basis vector along yarn direction
            !	e_yarn_2 			-> basis vector perpendicular to yarn direction
            !	
            !
            !	Outputs:-
            !   rotation_matrix		-> Rotation matrix
            !=======================================================================
            IMPLICIT NONE
            !Dummy variable declarations
            integer(kind=AbqIK),intent(in)::k,l
            integer(kind=AbqIK)::index_1
            
            !Local variable declarations
            
            
            !Variable initializations
            
            if (k <= l) then
				index_1 = k
            else 
				index_1 = k - l
            end if
                
     End Function
     
     Function rotate_fourth_order_tensor(original_tensor,rotation_matrix,NTENS) result(rotated_tensor)
            !=======================================================================
            !Function to compute the rotation matrix using yarn orientation.
            !
            !	Inputs:-
            !	e_yarn_1			-> basis vector along yarn direction
            !	e_yarn_2 			-> basis vector perpendicular to yarn direction
            !	
            !
            !	Outputs:-
            !   rotation_matrix		-> Rotation matrix
            !=======================================================================
             IMPLICIT NONE
            !Dummy variable declarations
            real(kind=AbqRK),intent(in)::original_tensor(:,:)
            real(kind=AbqRK),intent(in)::rotation_matrix(:,:)
            integer(kind=AbqIK),intent(in)::NTENS
            real(kind=AbqRK)::rotated_tensor(NTENS,NTENS)
            
            !Local variable declarations
            real(kind=AbqRK)::K_M(NTENS,NTENS), K_M_inv(NTENS,NTENS)	!Transformation matrix
            real(kind=AbqRK)::K_M_temp(NTENS) 	!Transformation matrix for strain tensor
            real(kind=AbqRK)::K_1(3,3),K_2(3,3),K_3(3,3),K_4(3,3) 	!Transformation matrix
            integer(kind=AbqIK)i,j !Running indicies
            
            !Variable initializations
            K_M		= 0.0_AbqRK
            K_M_inv = 0.0_AbqRK
            
            !K_1
            forall(i=1:3,j=1:3)K_1(i,j)= rotation_matrix(i,j)*rotation_matrix(i,j)
            
            !K_2
            forall(i=1:3,j=1:3)K_2(i,j)= rotation_matrix(i,mod_r(j+1,3))*&
										 rotation_matrix(i,mod_r(j+2,3))
            
            !K_3
            forall(i=1:3,j=1:3)K_3(i,j)= rotation_matrix(mod_r(i+1,3),j)*&
										 rotation_matrix(mod_r(i+2,3),j)
            
            !K_4
            forall(i=1:3,j=1:3)K_4(i,j)= rotation_matrix(mod_r(i+1,3),mod_r(j+1,3))*&
										 rotation_matrix(mod_r(i+2,3),mod_r(j+2,3)) +&
										 rotation_matrix(mod_r(i+1,3),mod_r(j+2,3))*&
										 rotation_matrix(mod_r(i+2,3),mod_r(j+1,3))
			!Assemble K_M
			K_M(:3,:3) 	= K_1(:,:)
			K_M(:3,4:6) 	= 2.0_AbqRK*K_2(:,:)
			K_M(4:6,:3) 	= K_3(:,:)
			K_M(4:6,4:6) 	= K_4(:,:)
			
			K_M_temp(:) = K_M(4,:)
            K_M(4,:) = K_M(6,:)
            K_M(6,:) = K_M_temp(:)
            
            K_M_temp(:) = K_M(:,4)
            K_M(:,4) = K_M(:,6)
            K_M(:,6) = K_M_temp(:)
            
            !compute inverse 
            call invert_matrix(K_M,K_M_inv)
			
			!rotate tensor
			rotated_tensor = matmul(matmul(K_M_inv,original_tensor),transpose(K_M_inv))
			
           
     End Function
     
	
	
END MODULE Yarn

!===============================================================================
!Module for matrix
!===============================================================================
MODULE Matrix
	USE ABQINTERFACE
    IMPLICIT NONE
    Private
    Public :: UMAT
    Contains
        !Include matrix material routine.
        INCLUDE 'UMATMatrix.f90'
        
END MODULE Matrix
