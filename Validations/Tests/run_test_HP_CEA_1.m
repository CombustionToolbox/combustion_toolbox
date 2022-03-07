function [max_rel_error_moles, max_rel_error_prop] = run_test_HP_CEA_1(value)
    % Run test_HP_CEA_1:
    % Contrasted with: NASA's Chemical Equilibrium with Applications software
    % Problem type: Adiabatic T and composition at constant p
    % Temperature [K]   = 300;
    % Pressure    [bar] = 1;
    % Equivalence ratio [-] = value
    % Initial mixture: Fuel + AIR_IDEAL (79% N2 + 21% O2)
    % List of species considered: ListSpecies('Soot Formation Extended')
    
    % Inputs
    Fuel = 'C2H2_acetylene';
    prefixDataName = 'C2H2';
    filename = {strcat(prefixDataName, '_air1_HP.out'), strcat(prefixDataName, '_air1_HP2.out'), strcat(prefixDataName, '_air1_HP3.out')};
    LS =  'Soot Formation Extended';
    DisplaySpecies = {'CO2', 'CO', 'H2O', 'H2', 'O2', 'N2', 'He', 'Ar',...
                      'HCN','H','OH','O','CN','NH3','CH4','C2H4','CH3',...
                      'NO','HCO','NH2','NH','N','CH','Cbgrb'};
    % Combustion Toolbox
    results_CT = run_CT('ListSpecies', LS, 'S_Fuel', Fuel,...
                        'S_Oxidizer', 'O2', 'S_Inert', 'N2',...
                        'EquivalenceRatio', value);
    % Load results CEA 
    results_CEA = data_CEA(filename, DisplaySpecies);
    % Compute error
    % * Molar fractions
    max_rel_error_moles = compute_error_moles_CEA(results_CT, results_CEA, 'phi', value, 'Xi', DisplaySpecies);
    % * Properties mixture 2
    max_rel_error_prop = compute_error_prop_CEA(results_CT, results_CEA, {'phi', 'phi', 'phi', 'phi', 'phi', 'phi', 'phi', 'phi'}, value, {'T', 'rho', 'h', 'e', 'g', 'cP', 'cV', 'gamma_s'}, 'mix2');
end