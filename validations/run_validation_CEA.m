function run_validation_CEA
    % Function that performs multiple tests to assess the quality and 
    % performance of the code. It compares the results against those
    % obtained from NASA's Chemical Equilibrium with Applications (CEA)
    % code

    % TP:
    run_validation_TP_CEA_1;
    run_validation_TP_CEA_2;
    run_validation_TP_CEA_3;
    run_validation_TP_CEA_4;
    run_validation_TP_CEA_6;

    % HP:
    run_validation_HP_CEA_1;
    run_validation_HP_CEA_2;
    run_validation_HP_CEA_3;
    run_validation_HP_CEA_4;

    % TV:
    run_validation_TV_CEA_1;
    run_validation_TV_CEA_2;

    % EV:
    run_validation_EV_CEA_1;

    % SHOCKS:
    run_validation_SHOCK_IONIZATION_CEA_1;
    run_validation_SHOCK_R_IONIZATION_CEA_1;

    % DET:
    run_validation_DET_CEA_1;
    run_validation_DET_CEA_2;
    run_validation_DET_CEA_3;
    run_validation_DET_CEA_4;
     
    % ROCKET:
    run_validation_ROCKET_CEA_17;
    run_validation_ROCKET_CEA_18;
    run_validation_ROCKET_CEA_20;
    run_validation_ROCKET_CEA_21;
    run_validation_ROCKET_CEA_22;
    run_validation_ROCKET_CEA_23;
    run_validation_ROCKET_CEA_24;
    run_validation_ROCKET_CEA_25;
end
