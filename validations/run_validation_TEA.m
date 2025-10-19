function run_validation_TEA
    % Function that performs multiple tests to assess the quality and
    % performance of the code. It compares the results against those 
    % obtained from the Thermochemical Equilibrium Abundances (TEA) code

    % Definitions
    FLAG_BENCHMARK = false;
    FLAG_EXPORT = true;

    % TP:
    run_validation_TP_TEA_1(FLAG_BENCHMARK, 'FLAG_EXPORT', FLAG_EXPORT);
    run_validation_TP_TEA_2(FLAG_BENCHMARK, 'FLAG_EXPORT', FLAG_EXPORT);
    run_validation_TP_TEA_3(FLAG_BENCHMARK, 'FLAG_EXPORT', FLAG_EXPORT);
    run_validation_TP_TEA_4(FLAG_BENCHMARK, 'FLAG_EXPORT', FLAG_EXPORT);
end