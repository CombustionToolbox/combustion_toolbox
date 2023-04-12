function [max_rel_error_moles, max_rel_error_prop] = run_test_SV_CEA_1(value, DB, DB_master)
    % Run test_TP_CEA_1:
    % Contrasted with: NASA's Chemical Equilibrium with Applications software
    % Problem type: Equilibrium composition at defined T and p
    % Density products [-] = value;
    % Pressure    [bar] = 10;
    % Equivalence ratio [-] = 0.5
    % Initial mixture: CH4 + AIR (78.084% N2, 20.9476% O2, 0.9365% Ar, 0.0319% CO2)
    % List of species considered: All (see routine find_products)
    
    % Inputs
    Fuel = {'CH4'};
    Oxidizer = {'N2', 'O2', 'Ar', 'CO2'};
    ratio_oxidizers_O2 = [78.084, 20.9476, 0.9365, 0.0319] ./ 20.9476;
    phi = 0.5;
    TR = 1500;
    pR = 10;
    prefixDataName = 'CH4_air';
    filename = {strcat(prefixDataName, '_SV1.out'), strcat(prefixDataName, '_SV2.out'), strcat(prefixDataName, '_SV3.out'), strcat(prefixDataName, '_SV4.out')};
    display_species = {'CO2', 'CO', 'H2O', 'H2', 'O2', 'N2',...
                       'HCN','H','OH','O','CN','NH3','CH4','C2H4','CH3',...
                       'NO','HCO','NH2','NH','N','CH','Cbgrb'};
    tolN = 1e-18;
    % Obtain volume ratio vP/vR
    mix1 = mixture(TR, pR, [Fuel, Oxidizer], [phi, 2 * ratio_oxidizers_O2], 'DB_master', DB_master, 'DB', DB);
    vP_vR = density(mix1) ./ value;

    % Combustion Toolbox
    results_CT = run_CT('ProblemType', 'SV',...
                        'TR', TR,...
                        'pR', pR,...
                        'species', 'Soot formation extended',...
                        'vP_vR', vP_vR,...
                        'S_Fuel', Fuel,...
                        'S_Oxidizer', Oxidizer,...
                        'ratio_oxidizers_O2', ratio_oxidizers_O2,...
                        'phi', phi,...
                        'tolN', tolN,...
                        'DB', DB,...
                        'DB_master', DB_master);
    % Load results CEA 
    results_CEA = data_CEA(filename, display_species);
    % Compute error
    % * Molar fractions
    max_rel_error_moles = compute_error_moles_CEA(results_CT, results_CEA, 'rho', value, 'Xi', display_species);
    % * Properties mixture 2
    properties_y = {'T', 'p', 'h', 'e', 'g', 'cP', 'cV', 'gamma_s', 'dVdT_p', 'dVdp_T', 'sound', 'W'};
    properties_x = create_cell_ntimes('rho', length(properties_y));
    max_rel_error_prop = compute_error_prop_CEA(results_CT, results_CEA, properties_x, value, properties_y, 'mix2');
end