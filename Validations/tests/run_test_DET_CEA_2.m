function [max_rel_error_moles, max_rel_error_prop] = run_test_DET_CEA_2(value)
    % Run test_DET_CEA_1:
    % Contrasted with: NASA's Chemical Equilibrium with Applications software
    % Problem type: Chapman-Jouguet Detonation
    % Temperature [K]   = 300;
    % Pressure    [bar] = 1;
    % Equivalence ratio [-] = 0.5:0.01:4
    % Initial mixture: Fuel + O2
    % List of species considered: ListSpecies('Soot Formation Extended')
    
    % Inputs
    Fuel = 'C2H2_acetylene';
    prefixDataName = 'C2H2';
    filename = {strcat(prefixDataName, '_O2_detonations.out'), strcat(prefixDataName, '_O2_detonations2.out'), strcat(prefixDataName, '_O2_detonations3.out')};

    LS =  'Soot Formation Extended';
    display_species = {'CO2', 'CO', 'H2O', 'H2', 'O2', 'N2', 'He', 'Ar',...
                      'HCN','H','OH','O','CN','NH3','CH4','C2H4','CH3',...
                      'NO','HCO','NH2','NH','N','CH','Cbgrb'};
    tolN = 1e-18;
    % Combustion Toolbox
    results_CT = run_CT('ProblemType', 'DET', ...
                        'pR', 1,...
                        'Species', LS, ...
                        'S_Fuel', Fuel, ...
                        'S_Oxidizer', {'O2'},...
                        'ratio_oxidizers_O2', 1,...
                        'EquivalenceRatio', value,...
                        'tolN', tolN);
    % Load results CEA 
    results_CEA = data_CEA(filename, display_species);
    % Compute error
    % * Molar fractions
    max_rel_error_moles = compute_error_moles_CEA(results_CT, results_CEA, 'phi', value, 'Xi', display_species);
    % * Properties mixture 2
    properties_y = {'T', 'p', 'rho', 'h', 'e', 'g', 'S', 'gamma_s', 'cP', 'cV', 'gamma_s', 'dVdT_p', 'dVdp_T', 'sound', 'W', 'u_preshock'};
    properties_x = {'phi', 'phi', 'phi', 'phi', 'phi', 'phi', 'phi', 'phi', 'phi', 'phi', 'phi', 'phi', 'phi', 'phi', 'phi', 'phi'};
    max_rel_error_prop = compute_error_prop_CEA(results_CT, results_CEA, properties_x, value, properties_y, 'mix2');
end