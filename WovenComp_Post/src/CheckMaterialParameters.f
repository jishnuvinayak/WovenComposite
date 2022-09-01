SUBROUTINE CheckMaterialParameters(PROPS)
    !===========================================================================
    !Dummy variable decelartions
    !===========================================================================
    Real(Kind=AbqRK),Intent(IN):: PROPS(:)
    !===========================================================================
    !Local variable declarations
    !===========================================================================
    Integer(Kind=AbqIK) :: material !Element material reference
    Integer(Kind=AbqIK) :: number_mat_parameters_input,number_mat_parameters_calculated ! Number of material parameters from input file , calculated
    Integer(Kind=AbqIK) :: number_kelvin_viogt_branches !Number of kelvin viogt branches used in the matrix material model
    !===========================================================================
    !Initialize variables
    !===========================================================================
    material=PROPS(1)
    number_mat_parameters_input = size(PROPS)-1 !-1 used to remove the material ref passed in along with parameters
    !===============================================================================
    !Subroutine logic -> Start
    !===============================================================================
    
    if(material == 2) then
        ! For yarn temp assumes 2. To be modfied after defining yarn material routine
        number_mat_parameters_calculated = 30
        if(number_mat_parameters_input /= number_mat_parameters_calculated) then
            call STDB_ABQERR(-3,"Error: Required number of material parameters missing for Yarn", 0, 0.0, " ")
        end if
    
    else if(material==1) then
        ! For matrix uses the relation : 2*N +19 from material model
        !number_kelvin_viogt_branches = PROPS(5)
        number_mat_parameters_calculated = 18 
        if(number_mat_parameters_input /= number_mat_parameters_calculated) then
            call STDB_ABQERR(-3,"Error: Required number of material parameters missing for Matrix", 0, 0.0, " ")
        end if
    else
        !Throw error message for undefined materials
        call STDB_ABQERR(-3,"Error: Material not defined", 0, 0.0, " ")
    end if

END SUBROUTINE
