function problems_solved = run_validation_HP_CEA_4
    % Run test validation_HP_CEA_4:
    % Contrasted with: NASA's Chemical Equilibrium with Applications software
    % Problem type: Adiabatic T and composition at constant p
    % Temperature [K]   = 300;
    % Pressure    [bar] = 1;
    % Equivalence ratio [-] = 0.5:0.01:4
    % Initial mixture: CH4 + O2
    % List of species considered:
    %  {'CO2', 'CO', 'H2O', 'H2', 'O2', 'N2', 'He', 'Ar',...
    %   'HCN','H','OH','O','CN','NH3','CH4','C2H4','CH3',...
    %   'NO','HCO','NH2','NH','N','CH','Cbgrb'}
    
    % Inputs
    Fuel = 'CH4';
    prefixDataName = Fuel;
    filename = {strcat(prefixDataName, '_O2_HP.out'), strcat(prefixDataName, '_O2_HP2.out'), strcat(prefixDataName, '_O2_HP3.out')};
    LS =  'Soot Formation Extended';
    DisplaySpecies = {'CO2', 'CO', 'H2O', 'H2', 'O2', 'H','OH','O',...
                      'CH4','C2H4','CH3','HCO','CH','Cbgrb'};
    tolN = 1e-18;
    % Combustion Toolbox
    results_CT = run_CT('ListSpecies', LS,...
                        'S_Fuel', Fuel,...
                        'S_Oxidizer', 'O2',...
                        'ratio_oxidizers_O2', 1,...
                        'EquivalenceRatio', 0.5:0.01:4,...
                        'tolN', tolN);
    problems_solved = length(results_CT.PD.range);
    % Load results CEA 
    results_CEA = data_CEA(filename, DisplaySpecies);
    % Display validation (plot)
    % * Molar fractions
    fig1 = plot_molar_fractions_validation(results_CT, results_CEA, 'phi', 'Xi', DisplaySpecies);
    % * Properties mixture 2
    fig2 = plot_properties_validation(results_CT, results_CEA, {'phi', 'phi', 'phi', 'phi', 'phi', 'phi', 'phi', 'phi'}, {'T', 'rho', 'h', 'e', 'g', 'cP', 'S', 'gamma_s'}, 'mix2');
    % Save plots
    folderpath = strcat(pwd,'\Validations\Figures\');
    stack_trace = dbstack;
    filename = stack_trace.name;
    saveas(fig1, strcat(folderpath, filename, '_molar'), 'svg');
    saveas(fig2, strcat(folderpath, filename, '_properties'), 'svg');
end