!Function 'GetNSTATEV' provides number of state variables needed for UMAT

FUNCTION GetNSTATV(NTENS,NDI,PROPS)
    !===========================================================================
    !Dummy variable decelartions
    !===========================================================================
    INTEGER(KIND=AbqIK):: GetNSTATV
    INTEGER(KIND=AbqIK), INTENT(IN):: NTENS,NDI
    REAL(KIND=AbqRK),INTENT(IN):: PROPS(:)
    !===========================================================================
    !Local variable decelartions
    !===========================================================================
    Integer (Kind=AbqIK)::material !Element material reference
    Integer (Kind=AbqIK)::number_kelvin_viogt_brances !Number of kelvin viogt branches used in the matrix material model
    !===========================================================================
    !Initialize variables
    !===========================================================================
    material=PROPS(1)
    
    !===============================================================================
    !Subroutine logic -> Start
    !===============================================================================
    if(material == 1) then
        !For matrix uses the relation : N*NTENS + NTENS+6 from material model
        number_kelvin_viogt_brances = PROPS(5)
        GetNSTATV = (NTENS*number_kelvin_viogt_brances) + (NTENS+6)
        GetNSTATV = 33
    else if(material == 2) then
        ! For yarn, temp returns 8.
        GetNSTATV = NTENS+2
        GetNSTATV = 8
    else
        !Throw error message for undefined materials
        call STDB_ABQERR(-3,"Error: Material not defined", 0, 0.0, " ")
    end if

END FUNCTION
