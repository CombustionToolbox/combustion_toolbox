function run_validation_TV_CEA_1
    % Run test validation_TV_CEA_1:
    % Contrasted with: NASA's Chemical Equilibrium with Applications software
    % Problem type: Equilibrium composition at defined T and v
    % Temperature [K]   = 3000;
    % Specific volume [m3/kg] = 0.8692;
    % Equivalence ratio [-] = 0.5:0.01:4
    % Initial mixture: CH4 + AIR (78.084% N2, 20.9476% O2, 0.9365% Ar, 0.0319% CO2)
    % List of species considered: list_species('Soot Formation Extended')

    % NOTE: VALIDATION HAS TO BE RECOMPUTED WITH A FIXED DENSITY VALUE!!!

    % Import packages
    import combustiontoolbox.databases.NasaDatabase
    import combustiontoolbox.core.*
    import combustiontoolbox.equilibrium.*
    
    % Benchmark?
    FLAG_BENCHMARK = false;

    % Definitions
    fuel = 'CH4';
    prefixDataName = fuel;
    for i = 36:-1:1
        filename{i} = strcat(prefixDataName, '_air_TV', sprintf('%d', i), '.out');
    end
    listSpecies =  'Soot Formation Extended';
    displaySpecies = {'CO2', 'CO', 'H2O', 'H2', 'O2', 'N2', 'Ar',...
                      'HCN','H','OH','O','CN','NH3','CH4','C2H4','CH3',...
                      'NO','HCO','NH2','NH','N','CH'};
    tolMoles = 1e-18;

    % Get Nasa database
    DB = NasaDatabase('FLAG_BENCHMARK', FLAG_BENCHMARK);
    
    % Define chemical system
    system = ChemicalSystem(DB, listSpecies);
    
    % Initialize mixture
    mix = Mixture(system);
    
    % Define chemical state
    set(mix, {fuel}, 'fuel', 1);
    set(mix, {'N2', 'O2', 'Ar', 'CO2'}, 'oxidizer', [78.084, 20.9476, 0.9365, 0.0319] / 20.9476);
    
    % Define properties
    mixArray = setProperties(mix, 'temperature', 3000, 'volume', 0.8692, 'equivalenceRatio', 0.5:0.01:4);
    
    % Initialize solver
    solver = EquilibriumSolver('problemType', 'TV', 'tolMoles', tolMoles, 'FLAG_RESULTS', false);
    
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
% fig2 = plot_properties_validation(results_CT, results_CEA, {'phi', 'phi', 'phi', 'phi', 'phi', 'phi', 'phi', 'phi'}, {'p', 'h', 'e', 'g', 'S', 'cP', 'cV', 'gamma_s'}, 'mix2');
% % Save plots
% folderpath = strcat(pwd,'\Validations\Figures\');
% stack_trace = dbstack;
% filename = stack_trace.name;
% saveas(fig1, strcat(folderpath, filename, '_molar'), 'svg');
% saveas(fig2, strcat(folderpath, filename, '_properties'), 'svg');