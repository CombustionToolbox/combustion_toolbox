function run_validation_TP_CEA_2
    % Run test validation_TP_CEA_2:
    % Contrasted with: NASA's Chemical Equilibrium with Applications software
    % Problem type: Equilibrium composition at defined T and p
    % Temperature [K]   = 2500;
    % Pressure    [bar] = 1;
    % Equivalence ratio [-] = 0.5:0.01:4
    % Initial mixture: C6H6 + O2
    % List of species considered: list_species('Soot Formation Extended')
    
    % Import packages
    import combustiontoolbox.databases.NasaDatabase
    import combustiontoolbox.core.*
    import combustiontoolbox.equilibrium.*
    import combustiontoolbox.utils.display.*
    
    % Benchmark?
    FLAG_BENCHMARK = false;

    % Definitions
    fuel = 'C6H6';
    prefixDataName = fuel;
    filename = {strcat(prefixDataName, '_O2_TP.out'), strcat(prefixDataName, '_O2_TP2.out'), strcat(prefixDataName, '_O2_TP3.out')};
    listSpecies = 'Soot Formation Extended';
    displaySpecies = {'CO2', 'CO', 'H2O', 'H2', 'O2', 'H','OH','O',...
                      'CH4','C2H4','CH3','HCO','CH','Cbgrb'};
    tolMoles = 1e-18;

    % Get Nasa database
    DB = NasaDatabase('FLAG_BENCHMARK', FLAG_BENCHMARK);
    
    % Define chemical system
    system = ChemicalSystem(DB, listSpecies);
    
    % Initialize mixture
    mix = Mixture(system);
    
    % Define chemical state
    set(mix, {fuel}, 'fuel', 1);
    set(mix, {'O2'}, 'oxidizer', 1);
    
    % Define properties
    mixArray = setProperties(mix, 'temperature', 2500, 'pressure', 1, 'equivalenceRatio', 0.5:0.01:4);
    
    % Initialize solver
    solver = EquilibriumSolver('problemType', 'TP', 'tolMoles', tolMoles, 'FLAG_RESULTS', false);
    
    % Solve problem
    solver.solveArray(mixArray);
    
    if FLAG_BENCHMARK
        return
    end

    % Load results CEA 
    resultsCEA = data_CEA(filename, displaySpecies);
    
    % Plot molar fractions
    plotComposition(mixArray(1), mixArray, 'equivalenceRatio', 'Xi', 'displaySpecies', displaySpecies, 'mintol', 1e-14, 'validation', resultsCEA);
end

% % Display validation (plot)
% % * Molar fractions
% [~, fig1] = plot_molar_fractions(results_CT, results_CT.PS.strP, 'phi', 'Xi', 'validation', results_CEA, 'display_species', display_species);
% % * Properties mixture 2
% fig2 = plot_properties_validation(results_CT, results_CEA, {'phi', 'phi', 'phi', 'phi', 'phi', 'phi', 'phi', 'phi'}, {'rho', 'h', 'e', 'g', 'S', 'cP', 'cV', 'gamma_s'}, 'mix2');
% % Save plots
% folderpath = strcat(pwd,'\Validations\Figures\');
% stack_trace = dbstack;
% filename = stack_trace.name;
% saveas(fig1, strcat(folderpath, filename, '_molar'), 'svg');
% saveas(fig2, strcat(folderpath, filename, '_properties'), 'svg');