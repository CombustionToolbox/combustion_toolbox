 function [max_rel_error_prop_mix1, max_rel_error_prop_mix2] = run_test_SHOCK_IONIZATION_CEA_1(value, DB, DB_master)
    % Run test_SHOCK_IONIZATION_CEA_1:
    % Contrasted with: NASA's Chemical Equilibrium with Applications software
    % Problem type: Planar incident shock wave
    % Temperature [K]   = 300
    % Pressure    [bar] = 1
    % Incident velocity [m/s] = value
    % Initial mixture: AIR (78.084% N2 + 20.9476% O2 + 0.9365% Ar + 0.0319% CO2)
    % List of species considered: list_species('Air_ions')
    
    % Inputs
    Oxidizer = {'N2', 'O2', 'Ar', 'CO2'};
    moles = [78.084, 20.9476, 0.9365, 0.0319] ./ 20.9476;
    LS = 'Air_ions';

    display_species = {'eminus','Ar','Arplus','C','Cplus','Cminus','CN','CNplus','CNminus',...
                      'CNN','CO','COplus','CO2','CO2plus',...
                      'N','Nplus','Nminus','NCO','NO','NOplus','NO2','NO2minus',...
                      'N2','N2plus','N2minus','NCN','N2O','N2Oplus','N3',...
                      'O','Oplus','Ominus','O2','O2plus','O2minus','O3'};
    filename = {'airNASA_incident_shocks_ionization1.out', 'airNASA_incident_shocks_ionization2.out',...
                'airNASA_incident_shocks_ionization3.out', 'airNASA_incident_shocks_ionization4.out',...
                'airNASA_incident_shocks_ionization5.out', 'airNASA_incident_shocks_ionization6.out'};
    %% TUNNING PARAMETERS
    tolN = 1e-14;
    % Combustion Toolbox
    results_CT = run_CT('ProblemType', 'SHOCK_I',...
                        'Species', LS,...
                        'S_Oxidizer', Oxidizer,...
                        'N_Oxidizer', moles,...
                        'u1', value,...
                        'tolN', tolN,...
                        'DB', DB,...
                        'DB_master', DB_master);
    % Load results CEA 
    results_CEA = data_CEA(filename, display_species, find_ind(list_species(LS), display_species));
    % Compute error
    % * Molar fractions
%     max_rel_error_moles = compute_error_moles_CEA(results_CT, results_CEA, 'u', value, 'Xi', DisplaySpecies);
    % * Properties mixture 1
    properties_y = {'T', 'p', 'h', 'e', 'g', 'S', 'cP', 'cV', 'gamma_s', 'sound', 'W'};
    properties_x = {'u', 'u', 'u', 'u', 'u', 'u', 'u', 'u', 'u', 'u', 'u'};
    max_rel_error_prop_mix1 = compute_error_prop_CEA(results_CT, results_CEA, properties_x, value, properties_y, 'mix1');
    % * Properties mixture 2
    properties_y = {'T', 'p', 'h', 'e', 'g', 'S', 'cP', 'cV', 'gamma_s', 'dVdT_p', 'dVdp_T', 'sound', 'W'};
    properties_x = {'u', 'u', 'u', 'u', 'u', 'u', 'u', 'u', 'u', 'u', 'u', 'u', 'u'};
    max_rel_error_prop_mix2 = compute_error_prop_CEA(results_CT, results_CEA, properties_x, value, properties_y, 'mix2');
end