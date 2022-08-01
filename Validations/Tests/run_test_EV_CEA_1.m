function [max_rel_error_moles, max_rel_error_prop] = run_test_EV_CEA_1(value)
    % Run test_EV_CEA_1:
    % Contrasted with: NASA's Chemical Equilibrium with Applications software
    % Problem type: Adiabatic T and composition at constant volume
    % Temperature [K]   = 300;
    % Pressure    [bar] = 1.01325;
    % Equivalence ratio [-] = value
    % Initial mixture: Fuel + AIR (78.084% N2, 20.9476% O2, 0.9365% Ar, 0.0319% CO2)
    % List of species considered: ListSpecies('Soot Formation Extended')
    
    % Inputs
    Fuel = 'CH4';
    prefixDataName = Fuel;
    for i = 36:-1:1
        filename{i} = strcat(prefixDataName, '_air_EV', sprintf('%d', i), '.out');
    end
    LS =  'Soot Formation Extended';
    DisplaySpecies = {'CO2', 'CO', 'H2O', 'H2', 'O2', 'N2', 'He', 'Ar',...
                      'HCN','H','OH','O','CN','NH3','CH4','C2H4','CH3',...
                      'NO','HCO','NH2','NH','N','CH','Cbgrb'};
    tolN = 1e-18;
    % Combustion Toolbox
    results_CT = run_CT('ProblemType', 'EV', ...
                        'pR', 1 * 1.01325,...
                        'Species', LS, ...
                        'S_Fuel', Fuel, ...
                        'S_Oxidizer', {'N2', 'O2', 'Ar', 'CO2'},...
                        'ratio_oxidizers_O2', [78.084, 20.9476, 0.9365, 0.0319] ./ 20.9476,...
                        'EquivalenceRatio', value,...
                        'tolN', tolN);
    % Load results CEA 
    results_CEA = data_CEA(filename, DisplaySpecies);
    % Compute error
    % * Molar fractions
    max_rel_error_moles = compute_error_moles_CEA(results_CT, results_CEA, 'phi', value, 'Xi', DisplaySpecies);
    % * Properties mixture 2
    properties_y = {'T', 'p', 'h', 'e', 'g', 'S', 'cP', 'cV', 'gamma_s', 'dVdT_p', 'dVdp_T', 'sound', 'W'};
    properties_x = {'phi', 'phi', 'phi', 'phi', 'phi', 'phi', 'phi', 'phi', 'phi', 'phi', 'phi', 'phi', 'phi'};
    max_rel_error_prop = compute_error_prop_CEA(results_CT, results_CEA, properties_x, value, properties_y, 'mix2');
end