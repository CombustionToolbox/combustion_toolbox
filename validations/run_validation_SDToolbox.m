function run_validation_SDToolbox
    % Function that performs multiple tests to assess the quality and
    % performance of the code. It compares the results against those
    % obtained from Caltech's Shock and Detonation Toolbox (SD-Toolbox)

    % Definitions
    FLAG_BENCHMARK = false;
    FLAG_EXPORT = true;

    % SHOCK_POLAR:
    run_validation_SHOCK_POLAR_SDToolbox_1(FLAG_BENCHMARK, 'FLAG_EXPORT', FLAG_EXPORT);
    run_validation_SHOCK_POLAR_SDToolbox_2(FLAG_BENCHMARK, 'FLAG_EXPORT', FLAG_EXPORT);

    % SHOCK_PRANDTL_MEYER
    run_validation_SHOCK_PRANDTL_MEYER_SDToolbox_1(FLAG_BENCHMARK, 'FLAG_EXPORT', FLAG_EXPORT);
    run_validation_SHOCK_PRANDTL_MEYER_SDToolbox_2(FLAG_BENCHMARK, 'FLAG_EXPORT', FLAG_EXPORT);
end