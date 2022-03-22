function problems_solved = run_validation_HP_CEA_1
    % Run test validation_HP_CEA_1:
    % Contrasted with: NASA's Chemical Equilibrium with Applications software
    % Problem type: Adiabatic T and composition at constant p
    % Temperature [K]   = 300;
    % Pressure    [bar] = 1;
    % Equivalence ratio [-] = 0.5:0.01:4
    % Initial mixture: C2H2_acetylene + AIR_IDEAL (79% N2 + 21% O2)
    % List of species considered: ListSpecies('Soot Formation Extended')
    
    % Inputs
    Fuel = 'C2H2_acetylene';
    prefixDataName = 'C2H2';
    filename = {strcat(prefixDataName, '_air1_HP.out'), strcat(prefixDataName, '_air1_HP2.out'), strcat(prefixDataName, '_air1_HP3.out')};
    LS =  'Soot Formation Extended';
    DisplaySpecies = {'CO2', 'CO', 'H2O', 'H2', 'O2', 'N2', 'He', 'Ar',...
                      'HCN','H','OH','O','CN','NH3','CH4','C2H4','CH3',...
                      'NO','HCO','NH2','NH','N','CH','Cbgrb'};
    tolN = 1e-18;
    % Combustion Toolbox
    results_CT = run_CT('ListSpecies', LS, 'S_Fuel', Fuel,...
                        'S_Oxidizer', 'O2', 'S_Inert', 'N2',...
                        'EquivalenceRatio', 0.5:0.01:4, 'tolN', tolN);
    problems_solved = length(results_CT.PD.range);
    % Load results CEA 
    results_CEA = data_CEA(filename, DisplaySpecies);
    % Display validation (plot)
    % * Molar fractions
    fig1 = plot_molar_fractions_validation(results_CT, results_CEA, 'phi', 'Xi', DisplaySpecies);
    % * Properties mixture 2
    fig2 = plot_properties_validation(results_CT, results_CEA, {'phi', 'phi', 'phi', 'phi', 'phi', 'phi', 'phi', 'phi'}, {'T', 'rho', 'h', 'e', 'g', 'cP', 'cV', 'gamma_s'}, 'mix2');
    % Save plots
    folderpath = strcat(pwd,'\Validations\Figures\');
    stack_trace = dbstack;
    filename = stack_trace.name;
    saveas(fig1, strcat(folderpath, filename, '_molar'), 'svg');
    saveas(fig2, strcat(folderpath, filename, '_properties'), 'svg');
end