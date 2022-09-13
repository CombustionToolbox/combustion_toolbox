function problems_solved = run_validation_DET_CEA_4
    % Run test validation_DET_CEA_4:
    % Contrasted with: NASA's Chemical Equilibrium with Applications software
    % Problem type: Chapman-Jouguet Detonation
    % Temperature [K]   = 300;
    % Pressure    [bar] = 1;
    % Equivalence ratio [-] = 0.5:0.01:4
    % Initial mixture: CH4 + O2
    % List of species considered: ListSpecies('Soot Formation Extended')
    
    % Inputs
    Fuel = 'CH4';
    prefixDataName = Fuel;
    filename = {strcat(prefixDataName, '_O2_detonations.out'), strcat(prefixDataName, '_O2_detonations2.out'), strcat(prefixDataName, '_O2_detonations3.out')};
    LS =  'Soot Formation Extended';
    display_species = {'CO2', 'CO', 'H2O', 'H2', 'O2', 'H','OH','O',...
                      'CH4','C2H4','CH3','HCO','CH'};
    tolN = 1e-18;
    % Combustion Toolbox
    results_CT = run_CT('ProblemType', 'DET',...
                        'ListSpecies', LS,...
                        'S_Fuel', Fuel,...
                        'S_Oxidizer', 'O2',...
                        'ratio_oxidizers_O2', 1,...
                        'EquivalenceRatio', 0.5:0.01:4,...
                        'tolN', tolN);
    problems_solved = length(results_CT.PD.range);
    % Load results CEA 
    results_CEA = data_CEA(filename, display_species);
    % Display validation (plot)
    % * Molar fractions
    fig1 = plot_molar_fractions_validation(results_CT, results_CEA, 'phi', 'Xi', display_species);
    % * Properties mixture 2 - 1
    fig2 = plot_properties_validation(results_CT, results_CEA, {'phi', 'phi', 'phi', 'phi', 'phi', 'phi', 'phi', 'phi'}, {'T', 'p', 'rho', 'h', 'e', 'g', 'S', 'gamma_s'}, 'mix2');
    % * Properties mixture 2 - 2
    fig3 = plot_properties_validation(results_CT, results_CEA, {'phi', 'phi', 'phi', 'phi', 'phi', 'phi', 'phi', 'phi'}, {'cP', 'cV', 'gamma_s', 'dVdT_p', 'dVdp_T', 'sound', 'W', 'u_preshock'}, 'mix2');
    % Save plots
    folderpath = strcat(pwd,'\Validations\Figures\');
    stack_trace = dbstack;
    filename = stack_trace.name;
    saveas(fig1, strcat(folderpath, filename, '_molar'), 'svg');
    saveas(fig2, strcat(folderpath, filename, '_properties_1'), 'svg');
    saveas(fig3, strcat(folderpath, filename, '_properties_2'), 'svg');
end