function problems_solved = run_validation_TV_CEA_1
    % Run test validation_TV_CEA_1:
    % Contrasted with: NASA's Chemical Equilibrium with Applications software
    % Problem type: Equilibrium composition at defined T and v
    % Temperature [K]   = 3000;
    % Pressure    [bar] = 1;
    % Equivalence ratio [-] = 0.5:0.01:4
    % Initial mixture: CH4 + AIR (78.084% N2, 20.9476% O2, 0.9365% Ar, 0.0319% CO2)
    % List of species considered: ListSpecies('Soot Formation Extended')
    
    % Inputs
    Fuel = 'CH4';
    prefixDataName = Fuel;
    for i = 36:-1:1
        filename{i} = strcat(prefixDataName, '_air_TV', sprintf('%d', i), '.out');
    end
    LS =  'Soot Formation Extended';
    DisplaySpecies = {'CO2', 'CO', 'H2O', 'H2', 'O2', 'N2', 'He', 'Ar',...
                      'HCN','H','OH','O','CN','NH3','CH4','C2H4','CH3',...
                      'NO','HCO','NH2','NH','N','CH','Cbgrb'};
    tolN = 1e-18;
    % Combustion Toolbox
    results_CT = run_CT('ProblemType', 'TV',...
                        'Temp', 3000,...
                        'pR', 1 * 1.01325,...
                        'Species', LS,...
                        'S_Fuel', Fuel,...
                        'S_Oxidizer', {'N2', 'O2', 'Ar', 'CO2'},...
                        'ratio_oxidizers_O2', [78.084, 20.9476, 0.9365, 0.0319] ./ 20.9476,...
                        'EquivalenceRatio', 0.5:0.01:4,...
                        'tolN', tolN);
    problems_solved = length(results_CT.PD.range);
    % Load results CEA 
    results_CEA = data_CEA(filename, DisplaySpecies);
    % Display validation (plot)
    % * Molar fractions
    fig1 = plot_molar_fractions_validation(results_CT, results_CEA, 'phi', 'Xi', DisplaySpecies);
    % * Properties mixture 2
    fig2 = plot_properties_validation(results_CT, results_CEA, {'phi', 'phi', 'phi', 'phi', 'phi', 'phi', 'phi', 'phi'}, {'p', 'h', 'e', 'g', 'S', 'cP', 'cV', 'gamma_s'}, 'mix2');
    % Save plots
    folderpath = strcat(pwd,'\Validations\Figures\');
    stack_trace = dbstack;
    filename = stack_trace.name;
    saveas(fig1, strcat(folderpath, filename, '_molar'), 'svg');
    saveas(fig2, strcat(folderpath, filename, '_properties'), 'svg');
end