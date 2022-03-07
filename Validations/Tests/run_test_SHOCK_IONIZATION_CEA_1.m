 function [max_rel_error_prop_mix1, max_rel_error_prop_mix2] = run_test_SHOCK_IONIZATION_CEA_1(value)
    % Run test_SHOCK_IONIZATION_CEA_1:
    % Contrasted with: NASA's Chemical Equilibrium with Applications software
    % Problem type: Planar incident shock wave
    % Temperature [K]   = 300
    % Pressure    [bar] = 1
    % Incident velocity [m/s] = value
    % Initial mixture: AIR (78.084% N2 + 20.9476% O2 + 0.9365% Ar + 0.0319% CO2)
    % List of species considered: ListSpecies('Air_ions')
    
    % Inputs
    Fuel = [];
    Oxidizer = 'O2';
    Inert    = {'N2', 'Ar', 'CO2'};
    proportion_inerts_O2 = [78.084, 0.9365, 0.0319] ./ 20.9476;
    LS =  {'eminus','Ar','Arplus','C','Cplus','Cminus','CN','CNplus','CNminus','CNN','CO','COplus','CO2','CO2plus','C2','C2plus','C2minus','CCN','CNC','OCCN','C2N2','C2O','C3','N','Nplus','Nminus','NCO','NO','NOplus','NO2','NO2minus','NO3','NO3minus','N2','N2plus','N2minus','N2O','NCN','N2Oplus','N2O3','N3','O','Oplus','Ominus','O2','O2plus','O2minus','O3'};

    DisplaySpecies = {'eminus','Ar','Arplus','C','Cplus','Cminus','CN','CNplus','CNminus',...
                      'CNN','CO','COplus','CO2','CO2plus',...
                      'N','Nplus','Nminus','NCO','NO','NOplus','NO2','NO2minus',...
                      'N2','N2plus','N2minus','NCN','N2O','N2Oplus','N3',...
                      'O','Oplus','Ominus','O2','O2plus','O2minus','O3'};
    filename = {'airNASA_incident_shocks_ionization1.out', 'airNASA_incident_shocks_ionization2.out',...
                'airNASA_incident_shocks_ionization3.out', 'airNASA_incident_shocks_ionization4.out',...
                'airNASA_incident_shocks_ionization5.out', 'airNASA_incident_shocks_ionization6.out'};
    %% TUNNING PARAMETERS
    tolN = 1e-16;
    % Combustion Toolbox
    results_CT = run_CT('ProblemType', 'SHOCK_I', 'Species', LS, 'S_Fuel', Fuel,...
                        'S_Oxidizer', Oxidizer, 'S_Inert', Inert,...
                        'proportion_inerts_O2', proportion_inerts_O2, 'u1', value,...
                        'tolN', tolN);
    % Load results CEA 
    results_CEA = data_CEA(filename, DisplaySpecies, find_ind(ListSpecies([], LS), DisplaySpecies));
    % Compute error
    % * Molar fractions
%     max_rel_error_moles = compute_error_moles_CEA(results_CT, results_CEA, 'u', value, 'Xi', DisplaySpecies);
    % * Properties mixture 1
    max_rel_error_prop_mix1 = compute_error_prop_CEA(results_CT, results_CEA, {'u', 'u', 'u', 'u', 'u', 'u', 'u', 'u'}, value, {'T', 'rho', 'h', 'e', 'g', 'cP', 'cV', 'gamma_s'}, 'mix1');
    % * Properties mixture 2
    max_rel_error_prop_mix2 = compute_error_prop_CEA(results_CT, results_CEA, {'u', 'u', 'u', 'u', 'u', 'u', 'u', 'u', 'u', 'u', 'u'}, value, {'T', 'p', 'rho', 'h', 'e', 'g', 'cP', 'gamma_s', 'cV', 'dVdT_p', 'dVdp_T'}, 'mix2');
end