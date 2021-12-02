function run_validation_TP_CEA_2
    % Run test validation_TP_CEA_2:
    % Contrasted with: NASA's Chemical Equilibrium with Applications software
    % Problem type: Equilibrium composition at defined T and p
    % Temperature [K]   = 2500;
    % Pressure    [bar] = 1;
    % Equivalence ratio [-] = 0.5:0.01:4
    % Initial mixture: C6H6 + O2
    % List of species considered: ListSpecies('Soot Formation Extended')
    
    % Inputs
    Fuel = 'C6H6';
    prefixDataName = Fuel;
    filename = {strcat(prefixDataName, '_O2_TP.out'), strcat(prefixDataName, '_O2_TP2.out'), strcat(prefixDataName, '_O2_TP3.out')};
    LS =  'Soot Formation Extended';
    DisplaySpecies = {'CO2', 'CO', 'H2O', 'H2', 'O2', 'N2', 'He', 'Ar',...
                      'HCN','H','OH','O','CN','NH3','CH4','C2H4','CH3',...
                      'NO','HCO','NH2','NH','N','CH','Cbgrb'};
    % Combustion Toolbox
    results_CT = run_CT('ProblemType', 'TP', 'Temp', 2500, 'Species', LS,...
                        'S_Fuel', Fuel,'S_Oxidizer', 'O2',...
                        'S_Inert', [], 'EquivalenceRatio', 0.5:0.01:4);
    % Load results CEA 
    results_CEA = data_CEA(filename, DisplaySpecies);
    % Display validation (plot)
    fig1 = plot_molar_fractions_validation(results_CT, results_CEA, 'phi', 'Xi', DisplaySpecies);
    % Save plots
    folderpath = strcat(pwd,'\Validations\Figures\');
    filename = 'validation_TP_CEA_2';
    saveas(fig1, strcat(folderpath, filename, '_molar'), 'svg');
end