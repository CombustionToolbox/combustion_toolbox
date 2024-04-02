function [max_rel_error_moles, max_rel_error_prop] = run_test_TP_CEA_1(value, DB)
    % Run test_TP_CEA_1:
    % Contrasted with: NASA's Chemical Equilibrium with Applications software
    % Problem type: Equilibrium composition at defined T and p
    % Temperature [K]   = value;
    % Pressure    [bar] = 1.01325;
    % Initial mixture: Si + 9 C6H5OH_phenol
    % List of species considered: All (see routine find_products)
    
    % Inputs
    Fuel = {'Si', 'C6H5OH_phenol'};
    moles_Fuel = [1, 9];
    prefixDataName = 'C6H5OH_phenol_and_Si';
    filename = {strcat(prefixDataName, '_TP1.out')};
    display_species = {'H2', 'H', 'CH4', 'C2H2_acetylene', 'SiC2',...
            'Si', 'C', 'C2H', 'C3', 'CO', 'CO2', 'Cbgrb', 'SiCbbb',...
            'SiO2ba_qzb','SiO2bb_qzb','SiO2bb_crtb','H2O','H2ObLb','H2Obcrb'};
    tolN = 1e-18;
    % Combustion Toolbox
    results_CT = run_CT('problemType', 'TP',...
                        'Temp', value,...
                        'Pressure', 1.01325,...
                        'listspecies', display_species,...
                        'S_Fuel', Fuel,...
                        'N_Fuel', moles_Fuel,...
                        'S_Oxidizer', [],...
                        'phi', [],...
                        'tolN', tolN,...
                        'DB', DB);
    % Load results CEA 
    results_CEA = data_CEA(filename, display_species);
    % Compute error
    % * Molar fractions
    max_rel_error_moles = compute_error_moles_CEA(results_CT, results_CEA, 'T', value, 'Xi', display_species);
    % * Properties mixture 2
    properties_y = {'p', 'h', 'e', 'g', 'S', 'cP', 'cV', 'gamma_s', 'dVdT_p', 'dVdp_T', 'sound', 'W'};
    properties_x = create_cell_ntimes('T', length(properties_y));
    max_rel_error_prop = compute_error_prop_CEA(results_CT, results_CEA, properties_x, value, properties_y, 'mix2');
end