function run_validation_Cantera
    % Function that performs multiple tests to assess the quality and
    % performance of the code. It compares the results against those 
    % obtained from the Cantera code

    % Definitions
    FLAG_BENCHMARK = false;
    FLAG_EXPORT = true;

    % TP:
    run_validation_TP_Cantera_1(FLAG_BENCHMARK, 'FLAG_EXPORT', FLAG_EXPORT);
end