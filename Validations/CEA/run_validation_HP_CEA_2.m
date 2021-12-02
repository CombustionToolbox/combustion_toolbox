function run_validation_HP_CEA_2
    % Run test validation_HP_CEA_2:
    % Contrasted with: NASA's Chemical Equilibrium with Applications software
    % Problem type: Adiabatic T and composition at constant p
    % Temperature [K]   = 300;
    % Pressure    [bar] = 1;
    % Equivalence ratio [-] = 0.5:0.01:4
    % Initial mixture: C2H2_acetylene + O2
    % List of species considered: ListSpecies('Soot Formation Extended')
    
    % Inputs
    Fuel = 'C2H2_acetylene';
    prefixDataName = 'C2H2';
    filename = {strcat(prefixDataName, '_O2_HP.out'), strcat(prefixDataName, '_O2_HP2.out'), strcat(prefixDataName, '_O2_HP3.out')};
    LS =  'Soot Formation Extended';
    DisplaySpecies = {'CO2', 'CO', 'H2O', 'H2', 'O2', 'N2', 'He', 'Ar',...
                      'HCN','H','OH','O','CN','NH3','CH4','C2H4','CH3',...
                      'NO','HCO','NH2','NH','N','CH','Cbgrb'};
    % Combustion Toolbox
    results_CT = run_CT('ListSpecies', LS, 'S_Fuel', Fuel,...
                        'S_Oxidizer', 'O2', 'S_Inert', [],...
                        'EquivalenceRatio', 0.5:0.01:4);
    % Load results CEA 
    results_CEA = data_CEA(filename, DisplaySpecies);
    % Display validation (plot)
    fig1 = plot_molar_fractions_validation(results_CT, results_CEA, 'phi', 'Xi', DisplaySpecies);
    % Save plots
    folderpath = strcat(pwd,'\Validations\Figures\');
    filename = 'validation_HP_CEA_2';
    saveas(fig1, strcat(folderpath, filename, '_molar'), 'svg');
end