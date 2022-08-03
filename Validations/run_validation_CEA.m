function run_validation_CEA
    % Function that runs several tests to evaluate the quality and
    % performance of the code compare with NASA's Chemical Equilibrium with
    % Applications software (CEA)
    
    % TP:
    run_validation_TP_CEA_1; 
    run_validation_TP_CEA_2;
    run_validation_TP_CEA_3;
    run_validation_TP_CEA_4;
    % HP:
    run_validation_HP_CEA_1; 
    run_validation_HP_CEA_2;
    run_validation_HP_CEA_3;
    run_validation_HP_CEA_4;
    % TV:
    run_validation_TV_CEA_1;
    % EV:
    run_validation_EV_CEA_1;
    % DET: 
    run_validation_DET_CEA_1;
    run_validation_DET_CEA_2;
    run_validation_DET_CEA_3;
    run_validation_DET_CEA_4;
    % SHOCKS: 
    run_validation_SHOCK_IONIZATION_CEA_1;
    run_validation_SHOCK_R_IONIZATION_CEA_1;
    % ROCKET:
    run_validation_ROCKET_CEA_17;
    run_validation_ROCKET_CEA_18;
end