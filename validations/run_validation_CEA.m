function run_validation_CEA
    % Function that performs multiple tests to assess the quality and 
    % performance of the code. It compares the results against those
    % obtained from NASA's Chemical Equilibrium with Applications (CEA)
    % code
    
    % Definitions
    FLAG_BENCHMARK = false;
    FLAG_EXPORT = true;

    % TP:
    run_validation_TP_CEA_1(FLAG_BENCHMARK, 'FLAG_EXPORT', FLAG_EXPORT);
    run_validation_TP_CEA_2(FLAG_BENCHMARK, 'FLAG_EXPORT', FLAG_EXPORT);
    run_validation_TP_CEA_3(FLAG_BENCHMARK, 'FLAG_EXPORT', FLAG_EXPORT);
    run_validation_TP_CEA_4(FLAG_BENCHMARK, 'FLAG_EXPORT', FLAG_EXPORT);
    run_validation_TP_CEA_6(FLAG_BENCHMARK, 'FLAG_EXPORT', FLAG_EXPORT);

    % HP:
    run_validation_HP_CEA_1(FLAG_BENCHMARK, 'FLAG_EXPORT', FLAG_EXPORT);
    run_validation_HP_CEA_2(FLAG_BENCHMARK, 'FLAG_EXPORT', FLAG_EXPORT);
    run_validation_HP_CEA_3(FLAG_BENCHMARK, 'FLAG_EXPORT', FLAG_EXPORT);
    run_validation_HP_CEA_4(FLAG_BENCHMARK, 'FLAG_EXPORT', FLAG_EXPORT);
    run_validation_HP_CEA_5(FLAG_BENCHMARK, 'FLAG_EXPORT', FLAG_EXPORT);

    % TV:
    run_validation_TV_CEA_1(FLAG_BENCHMARK, 'FLAG_EXPORT', FLAG_EXPORT);
    run_validation_TV_CEA_2(FLAG_BENCHMARK, 'FLAG_EXPORT', FLAG_EXPORT);

    % EV:
    run_validation_EV_CEA_1(FLAG_BENCHMARK, 'FLAG_EXPORT', FLAG_EXPORT);

    % SHOCKS:
    run_validation_SHOCK_IONIZATION_CEA_1(FLAG_BENCHMARK, 'FLAG_EXPORT', FLAG_EXPORT);
    run_validation_SHOCK_R_IONIZATION_CEA_1(FLAG_BENCHMARK, 'FLAG_EXPORT', FLAG_EXPORT);

    % DET:
    run_validation_DET_CEA_1(FLAG_BENCHMARK, 'FLAG_EXPORT', FLAG_EXPORT);
    run_validation_DET_CEA_2(FLAG_BENCHMARK, 'FLAG_EXPORT', FLAG_EXPORT);
    run_validation_DET_CEA_3(FLAG_BENCHMARK, 'FLAG_EXPORT', FLAG_EXPORT);
    run_validation_DET_CEA_4(FLAG_BENCHMARK, 'FLAG_EXPORT', FLAG_EXPORT);
     
    % ROCKET:
    run_validation_ROCKET_CEA_17(FLAG_BENCHMARK, 'FLAG_EXPORT', FLAG_EXPORT);
    run_validation_ROCKET_CEA_18(FLAG_BENCHMARK, 'FLAG_EXPORT', FLAG_EXPORT);
    run_validation_ROCKET_CEA_20(FLAG_BENCHMARK, 'FLAG_EXPORT', FLAG_EXPORT);
    run_validation_ROCKET_CEA_21(FLAG_BENCHMARK, 'FLAG_EXPORT', FLAG_EXPORT);
    run_validation_ROCKET_CEA_22(FLAG_BENCHMARK, 'FLAG_EXPORT', FLAG_EXPORT);
    run_validation_ROCKET_CEA_23(FLAG_BENCHMARK, 'FLAG_EXPORT', FLAG_EXPORT);
    run_validation_ROCKET_CEA_24(FLAG_BENCHMARK, 'FLAG_EXPORT', FLAG_EXPORT);
    run_validation_ROCKET_CEA_25(FLAG_BENCHMARK, 'FLAG_EXPORT', FLAG_EXPORT);
end
