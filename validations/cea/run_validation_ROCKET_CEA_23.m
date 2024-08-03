function run_validation_ROCKET_CEA_23
    % Run test validation_ROCKET_CEA_23:
    % Contrasted with: NASA's Chemical Equilibrium with Applications software
    % Problem type: Equilibrium composition at exit of the rocket nozzle
    % Pressure chamber [bar] = 101.325;
    % Model: Finite Area Chamber (FAC)
    % Area ratio A_c/A_t = 2;
    % Area ratio A_e/A_t = 3;
    % Equivalence ratio [-] = 0.5:0.01:4
    % Initial mixture: CH6N2bLb + MMH
    % List of species considered: list_species('Soot formation extended')

    % Import packages
    import combustiontoolbox.databases.NasaDatabase
    import combustiontoolbox.core.*
    import combustiontoolbox.rocket.*
    import combustiontoolbox.utils.display.*
    
    % Benchmark?
    FLAG_BENCHMARK = false;

    % Definitions
    fuel = 'CH6N2bLb';
    areaRatioChamber = 2;
    areaRatio = 3;
    prefixDataName = 'MMH';
    filename = {strcat(prefixDataName, '_N2O4_ROCKET_FAC1.out'), strcat(prefixDataName, '_N2O4_ROCKET_FAC2.out')};
    listSpecies = 'Soot formation extended';
    displaySpecies = {'N2', 'H2O', 'O2', 'CO2', 'NO', 'OH', 'O', 'CO',...
        'H2', 'NO2', 'HO2', 'H', 'H2O2', 'N2O', 'HNO'};
    tolMoles = 1e-18;

    % Get Nasa database
    DB = NasaDatabase('FLAG_BENCHMARK', FLAG_BENCHMARK);
    
    % Define chemical system
    system = ChemicalSystem(DB, listSpecies);
    
    % Initialize mixture
    mix = Mixture(system);
    
    % Define chemical state
    set(mix, {fuel}, 'fuel', 1);
    set(mix, {'N2O4'}, 'oxidizer', 1);
    
    % Define properties
    mixArray1 = setProperties(mix, 'temperature', 300, 'pressure', 101.325, 'equivalenceRatio', 0.5:0.01:4, 'areaRatioChamber', areaRatioChamber, 'areaRatio', areaRatio);
    
    % Initialize solver
    solver = RocketSolver('problemType', 'ROCKET_FAC', 'tolMoles', tolMoles, 'FLAG_RESULTS', false);
    
    % Solve problem
    [~, ~, ~, ~, mixArray4] = solver.solveArray(mixArray1);

    if FLAG_BENCHMARK
        return
    end
    
    % Load results CEA 
    resultsCEA = data_CEA(filename, displaySpecies);
    
    % Plot molar fractions
    fig1 = plotComposition(mixArray4(1), mixArray4, 'equivalenceRatio', 'Xi', 'displaySpecies', displaySpecies, 'mintol', 1e-14, 'validation', resultsCEA);

    % Plot properties
    fig2 = plotProperties(repmat({'equivalenceRatio'}, 1, 12), mixArray4, {'T', 'rho', 'h', 'e', 'g', 'cp', 's', 'gamma_s', 'sound', 'u', 'I_sp', 'I_vac'}, mixArray4, 'basis', {[], [], 'mi', 'mi', 'mi', 'mi', 'mi', [], [], [], [], []}, 'validation', resultsCEA);

    % Save plots
    folderpath = fullfile(pwd, 'validations', 'figures');
    stack_trace = dbstack;
    filename = stack_trace.name;
    saveas(fig1, fullfile(folderpath, strcat(filename, '_molar')), 'svg');
    saveas(fig2, fullfile(folderpath, strcat(filename, '_properties')), 'svg');
end