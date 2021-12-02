function run_validations_CEA
    % Function that runs several tests to evaluate the quality and
    % performance of the code compare with NASA's Chemical Equilibrium with
    % Applications software (CEA)
    
    % TP:
    run_validation_TP_1; 
    run_validation_TP_2;
    run_validation_TP_3;
    run_validation_TP_4;
    % HP:
    run_validation_HP_1; 
    run_validation_HP_2;
    run_validation_HP_3;
    run_validation_HP_4;
    % DET: 
    run_validation_DET_1;
    run_validation_DET_2;
    run_validation_DET_3;
    run_validation_DET_4;
end