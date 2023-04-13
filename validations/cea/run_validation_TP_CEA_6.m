function problems_solved = run_validation_TP_CEA_6
    % Run test validation_TP_CEA_6:
    % Contrasted with: NASA's Chemical Equilibrium with Applications software
    % Problem type: Equilibrium composition at defined T and p
    % Temperature [K]   = 2000:50:5000;
    % Pressure    [bar] = 1.01325;
    % Initial mixture: Si + 9 C6H5OH_phenol
    % List of species considered: All (see routine find_products)
    %
    % This example is obtained from [1]
    % 
    % [1] Scoggins, J. B., & Magin, T. E. (2015). Gibbs function continuation
    %     for linearly constrained multiphase equilibria. Combustion and Flame,
    %     162(12), 4514-4522.
    
    % Inputs
    Fuel = {'Si', 'C6H5OH_phenol'};
    moles_Fuel = [1, 9];
    prefixDataName = 'C6H5OH_phenol_and_Si';
    filename = {strcat(prefixDataName, '_TP1.out')};
    display_species = {'H2', 'H', 'CH4', 'C2H2_acetylene', 'SiC2',...
            'Si', 'C', 'C2H', 'C3', 'CO', 'CO2', 'Cbgrb', 'SiCbbb',...
            'SiO2ba_qzb','SiO2bb_qzb','SiO2bb_crtb','H2O','H2ObLb','H2Obcrb'};
    tolN = 1e-18;
    % Load results CEA 
    results_CEA = data_CEA(filename, display_species);
    % Combustion Toolbox
    results_CT = run_CT('ProblemType', 'TP',...
                        'Temp', 200:10:5000,...
                        'Pressure', 1.01325,...
                        'S_Fuel', Fuel,...
                        'N_Fuel', moles_Fuel,...
                        'S_Oxidizer', [],...
                        'phi', [],...
                        'tolN', tolN);
    problems_solved = length(results_CT.PD.range);
    % Display validation (plot)
    % * Molar fractions
    [~, fig1] = plot_molar_fractions(results_CT, results_CT.PS.strP, 'T', 'Xi', 'validation', results_CEA, 'display_species', display_species);
    % * Properties mixture 2
    fig2 = plot_properties_validation(results_CT, results_CEA, {'T', 'T', 'T', 'T', 'T', 'T', 'T', 'T', 'T', 'T'}, {'rho', 'h', 'e', 'g', 'S', 'cP', 'cV', 'gamma_s', 'dVdp_T', 'dVdT_p'}, 'mix2');
    % Save plots
    folderpath = strcat(pwd,'\Validations\Figures\');
    stack_trace = dbstack;
    filename = stack_trace.name;
    saveas(fig1, strcat(folderpath, filename, '_molar'), 'svg');
    saveas(fig2, strcat(folderpath, filename, '_properties'), 'svg');
end