function run_validation_TV_CEA_2
    % Run test validation_TV_CEA_2:
    % Contrasted with: NASA's Chemical Equilibrium with Applications software
    % Problem type: Equilibrium composition at defined T and v
    % Temperature [K]   = 1200;
    % Specific volume [m3/kg] = 1;
    % Equivalence ratio [-] = 0.5:0.01:4
    % Initial mixture: Jet_AbLb + AIR (78.084% N2, 20.9476% O2, 0.9365% Ar, 0.0319% CO2)
    % List of species considered: All (see method findProducts from ChemicalSystem class)
    
    % Import packages
    import combustiontoolbox.databases.NasaDatabase
    import combustiontoolbox.core.*
    import combustiontoolbox.equilibrium.*
    import combustiontoolbox.utils.display.*
    
    % Benchmark?
    FLAG_BENCHMARK = false;

    % Definitions
    fuel = 'Jet_AbLb';
    prefixDataName = fuel;
    
    for i = 2:-1:1
        filename{i} = strcat(prefixDataName, '_air_TV', sprintf('%d', i), '.out');
    end

    displaySpecies = {'CO2', 'CO', 'H2O', 'H2', 'O2', 'N2', 'Ar',...
                      'HCN','H','OH','O','CN','NH3','CH4','C2H4', 'CH3',...
                      'NO', 'HCO', 'NH2', 'NH','N', 'CH', 'Cbgrb'};
    tolMoles = 1e-30;

    % Get Nasa database
    DB = NasaDatabase('FLAG_BENCHMARK', FLAG_BENCHMARK);
    
    % Define chemical system
    system = ChemicalSystem(DB);
    
    % Initialize mixture
    mix = Mixture(system);
    
    % Define chemical state
    set(mix, {fuel}, 'fuel', 1);
    set(mix, {'N2', 'O2', 'Ar', 'CO2'}, 'oxidizer', [78.084, 20.9476, 0.9365, 0.0319] / 20.9476);
    
    % Define properties
    mixArray = setProperties(mix, 'temperature', 1200, 'volume', 1, 'equivalenceRatio', 0.5:0.01:4);
    
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
    fig1 = plotComposition(mixArray(1), mixArray, 'equivalenceRatio', 'Xi', 'displaySpecies', displaySpecies, 'mintol', 1e-14, 'validation', resultsCEA);

    % Plot properties
    fig2 = plotProperties(repmat({'equivalenceRatio'}, 1, 10), mixArray, {'p', 'h', 'e', 'g', 'cp', 's', 'gamma_s', 'sound'}, mixArray, 'basis', {[], 'mi', 'mi', 'mi', 'mi', 'mi', [], []}, 'validation', resultsCEA);

    % Save plots
    folderpath = fullfile(pwd, 'validations', 'figures');
    stack_trace = dbstack;
    filename = stack_trace.name;
    saveas(fig1, fullfile(folderpath, strcat(filename, '_molar')), 'svg');
    saveas(fig2, fullfile(folderpath, strcat(filename, '_properties')), 'svg');
end