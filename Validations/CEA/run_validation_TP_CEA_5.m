function problems_solved = run_validation_TP_CEA_5
    % Run test validation_TP_CEA_5:
    % Contrasted with: NASA's Chemical Equilibrium with Applications software
    % Problem type: Equilibrium composition at defined T and p
    % Temperature [K]   = 2000:100:20000;
    % Pressure    [bar] = 1;
    % Equivalence ratio [-] = 1.1
    % Initial mixture: C6H6 + AIR_IDEAL (79% N2 + 21% O2)
    % List of species considered: ListSpecies('Soot Formation Extended')
    
    % Inputs
    Fuel = 'C6H6';
    prefixDataName = Fuel;
    filename = {strcat(prefixDataName, '_air_T_TP1.out'), strcat(prefixDataName, '_air_T_TP2.out')};
    LS =  'Soot Formation Extended';
    DisplaySpecies = {'CO2', 'CO', 'H2O', 'H2', 'O2', 'N2', 'He', 'Ar',...
                      'HCN','H','OH','O','CN','NH3','CH4','C2H4','CH3',...
                      'NO','HCO','NH2','NH','N','CH','Cbgrb'};
    tolN = 1e-18;
    % Combustion Toolbox
%     results_CT = run_CT('ProblemType', 'TP', 'Temp', 2000:100:20000, 'Species', LS,...
%                         'S_Fuel', Fuel,'S_Oxidizer', 'O2',...
%                         'S_Inert', 'N2', 'EquivalenceRatio', 1.1,...
%                         'tolN', tolN);
    problems_solved = length(results_CT.PD.range);
    % Load results CEA 
    results_CEA = data_CEA(filename, DisplaySpecies);
    % Display validation (plot)
    % * Molar fractions
    fig1 = plot_molar_fractions_validation(results_CT, results_CEA, 'phi', 'Xi', DisplaySpecies);
    % * Properties mixture 2
    fig2 = plot_properties_validation(results_CT, results_CEA, {'phi', 'phi', 'phi', 'phi', 'phi', 'phi', 'phi', 'phi'}, {'rho', 'h', 'e', 'g', 'S', 'cP', 'cV', 'gamma_s'}, 'mix2');
    % Save plots
    folderpath = strcat(pwd,'\Validations\Figures\');
    stack_trace = dbstack;
    filename = stack_trace.name;
    saveas(fig1, strcat(folderpath, filename, '_molar'), 'svg');
    saveas(fig2, strcat(folderpath, filename, '_properties'), 'svg');
end