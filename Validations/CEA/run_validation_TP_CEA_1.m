function problems_solved = run_validation_TP_CEA_1
    % Run test validation_TP_CEA_1:
    % Contrasted with: NASA's Chemical Equilibrium with Applications software
    % Problem type: Equilibrium composition at defined T and p
    % Temperature [K]   = 2500;
    % Pressure    [bar] = 1;
    % Equivalence ratio [-] = 0.5:0.01:4
    % Initial mixture: C6H6 + AIR_IDEAL (79% N2 + 21% O2)
    % List of species considered: ListSpecies('Soot Formation Extended')
    
    % Inputs
    Fuel = 'C6H6';
    prefixDataName = Fuel;
    filename = {strcat(prefixDataName, '_air1_TP.out'), strcat(prefixDataName, '_air1_TP2.out'), strcat(prefixDataName, '_air1_TP3.out')};
    LS = 'Soot Formation Extended';
    display_species = {'CO2', 'CO', 'H2O', 'H2', 'O2', 'N2',...
                      'HCN','H','OH','O','CN','NH3','CH4','C2H4','CH3',...
                      'NO','HCO','NH2','NH','N','CH','Cbgrb'};
    tolN = 1e-18;
    % Combustion Toolbox
    results_CT = run_CT('ProblemType', 'TP',...
                        'Temp', 2500,...
                        'Species', LS,...
                        'S_Fuel', Fuel,...
                        'S_Oxidizer', {'N2', 'O2'},...
                        'ratio_oxidizer_O2', [79, 21]/21,...
                        'EquivalenceRatio', 0.5:0.01:4,...
                        'tolN', tolN);
    problems_solved = length(results_CT.PD.range);
    % Load results CEA 
    results_CEA = data_CEA(filename, display_species);
    % Display validation (plot)
    % * Molar fractions
    fig1 = plot_molar_fractions_validation(results_CT, results_CEA, 'phi', 'Xi', display_species);
    % * Properties mixture 2
    fig2 = plot_properties_validation(results_CT, results_CEA, {'phi', 'phi', 'phi', 'phi', 'phi', 'phi', 'phi', 'phi'}, {'rho', 'h', 'e', 'g', 'S', 'cP', 'cV', 'gamma_s'}, 'mix2');
    % Save plots
    folderpath = strcat(pwd,'\Validations\Figures\');
    stack_trace = dbstack;
    filename = stack_trace.name;
    saveas(fig1, strcat(folderpath, filename, '_molar'), 'svg');
    saveas(fig2, strcat(folderpath, filename, '_properties'), 'svg');
end