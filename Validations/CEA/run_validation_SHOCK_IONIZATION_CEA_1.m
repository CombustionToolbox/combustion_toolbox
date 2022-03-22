 function run_validation_SHOCK_IONIZATION_CEA_1
    % Run test validation_SHOCK_IONIZATION_CEA_1:
    % Contrasted with: NASA's Chemical Equilibrium with Applications software
    % Problem type: Planar incident shock wave
    % Temperature [K]   = 300
    % Pressure    [bar] = 1
    % Incident velocity [m/s] = [~308, 15000]
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
    u1 = logspace(2, 5, 500); u1 = u1(u1<13000); u1 = u1(u1>=357);
    filename = {'airNASA_incident_shocks_ionization1.out', 'airNASA_incident_shocks_ionization2.out',...
                'airNASA_incident_shocks_ionization3.out', 'airNASA_incident_shocks_ionization4.out',...
                'airNASA_incident_shocks_ionization5.out', 'airNASA_incident_shocks_ionization6.out'};
    % TUNNING PARAMETERS
    tolN = 1e-14;
    % Combustion Toolbox
    results_CT = run_CT('ProblemType', 'SHOCK_I', 'Species', LS, 'S_Fuel', Fuel,...
                        'S_Oxidizer', Oxidizer, 'S_Inert', Inert,...
                        'proportion_inerts_O2', proportion_inerts_O2, 'u1', u1,...
                        'tolN', tolN);
    % Load results CEA 
    results_CEA = data_CEA(filename, DisplaySpecies, find_ind(ListSpecies([], LS), DisplaySpecies));
    % Display validations (plot)
    % * Molar fractions
%     fig1 = plot_molar_fractions_validation(results_CT, results_CEA.mix2, 'T', 'Xi', DisplaySpecies);
    % * Properties mixture 1
    fig2 = plot_properties_validation(results_CT, results_CEA.mix1, {'u', 'u', 'u', 'u', 'u', 'u', 'u', 'u'}, {'T', 'p', 'rho', 'h', 'e', 'g', 'cP', 'gamma_s'}, 'mix1');
    % * Properties mixture 2 - 1
    fig3 = plot_properties_validation(results_CT, results_CEA.mix2, {'u', 'u', 'u', 'u', 'u', 'u', 'u', 'u'}, {'T', 'p', 'rho', 'h', 'e', 'g', 'cP', 'gamma_s'}, 'mix2');
    % * Properties mixture 2 - 2
    fig4 = plot_properties_validation(results_CT, results_CEA.mix2, {'u', 'u', 'u'}, {'cV', 'dVdT_p', 'dVdp_T'}, 'mix2');
    % Save plots
    folderpath = strcat(pwd,'\Validations\Figures\');
    stack_trace = dbstack;
    filename = stack_trace.name;
%     saveas(fig1, strcat(folderpath, filename, '_molar'), 'svg');
    saveas(fig2, strcat(folderpath, filename, '_properties_1'), 'svg');
    saveas(fig3, strcat(folderpath, filename, '_properties_2'), 'svg');
    saveas(fig4, strcat(folderpath, filename, '_properties_3'), 'svg');
end