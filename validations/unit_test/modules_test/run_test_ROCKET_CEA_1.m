function [max_rel_error_moles, max_rel_error_prop] = run_test_ROCKET_CEA_1(value, DB, DB_master)
    % Run test run_test_ROCKET_CEA_1:
    % Contrasted with: NASA's Chemical Equilibrium with Applications software
    % Problem type: Equilibrium composition at exit of the rocket nozzle
    % Pressure chamber [bar] = 101.325;
    % Model: Finite Area Chamber (FAC)
    % Area ratio A_c/A_t = 2;
    % Area ratio A_e/A_t = 3;
    % Equivalence ratio [-] = 0.5:0.01:4
    % Initial mixture: LCH4 + LOX === CH4bLb + O2bLb
    % List of species considered: list_species('Soot Formation Extended')
    
    % Inputs
    Fuel = 'CH4bLb';
    Aratio_c = 2;
    Aratio = 3;
    prefixDataName = 'LCH4';
    filename = {[prefixDataName, '_LOX_ROCKET_FAC1.out'], [prefixDataName, '_LOX_ROCKET_FAC2.out']};
    LS =  'Soot Formation Extended';
    display_species = {'CO2','CO','H2O','H2','O2','C2H2_acetylene',...
          'C2H4','C2H6','CH2CO_ketene','CH3','CH3CHO_ethanal','CH3OH',...
          'CH4','COOH','H','H2O2','HCHO_formaldehy','HCO','HCOOH','HO2',...
          'O','OH','Cbgrb'};
    tolN = 1e-18;
    % Combustion Toolbox
    results_CT = run_CT('ProblemType', 'ROCKET',...
                        'TR', 298.15,...
                        'pR', 101.325,...
                        'Species', LS,...
                        'S_Fuel', Fuel,...
                        'S_Oxidizer', 'O2bLb',...
                        'ratio_oxidizers_O2', 1,...
                        'EquivalenceRatio', value,...
                        'tolN', tolN,...
                        'FLAG_IAC', false,...
                        'Aratio_c', Aratio_c,...
                        'Aratio', Aratio,...
                        'DB', DB,...
                        'DB_master', DB_master);
    % Load results CEA 
    results_CEA = data_CEA(filename, display_species);
    % Compute error
    % * Molar fractions
    max_rel_error_moles = compute_error_moles_CEA(results_CT, results_CEA, 'phi', value, 'Xi', display_species);
    % * Properties mixture 2
    properties_y = {'T', 'p', 'h', 'e', 'g', 'S', 'cP', 'cV', 'gamma_s', 'dVdT_p', 'dVdp_T', 'sound', 'W', 'u', 'I_sp', 'I_vac'};
    properties_x = create_cell_ntimes('phi', length(properties_y));
    max_rel_error_prop = compute_error_prop_CEA(results_CT, results_CEA, properties_x, value, properties_y, 'mix2');
end