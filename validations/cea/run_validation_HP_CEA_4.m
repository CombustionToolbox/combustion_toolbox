function metadata = run_validation_HP_CEA_4(varargin)
    % Run test validation_HP_CEA_4:
    % Contrasted with: NASA's Chemical Equilibrium with Applications software
    % Problem type: Adiabatic T and composition at constant p
    % Temperature [K]   = 300;
    % Pressure    [bar] = 1;
    % Equivalence ratio [-] = 0.5:0.01:4
    % Initial mixture: CH4 + O2
    % List of species considered: list_species('Soot Formation Extended')

    % Import packages
    import combustiontoolbox.databases.NasaDatabase
    import combustiontoolbox.core.*
    import combustiontoolbox.equilibrium.*
    import combustiontoolbox.utils.display.*
    
    % Default values
    metadata = [];
    DEFAULT_FLAG_BENCHMARK = false;
    DEFAULT_FLAG_EXPORT = false;

    % Input parser
    p = inputParser;
    addOptional(p, 'FLAG_BENCHMARK', DEFAULT_FLAG_BENCHMARK, @(x) islogical(x));
    addParameter(p, 'FLAG_EXPORT', DEFAULT_FLAG_EXPORT, @(x) islogical(x));
    parse(p, varargin{:});
    FLAG_BENCHMARK = p.Results.FLAG_BENCHMARK;
    FLAG_EXPORT = p.Results.FLAG_EXPORT;

    % Definitions
    fuel = 'CH4';
    prefixDataName = fuel;
    filename = {strcat(prefixDataName, '_O2_HP.out'), strcat(prefixDataName, '_O2_HP2.out'), strcat(prefixDataName, '_O2_HP3.out')};
    listSpecies =  'Soot Formation Extended';
    displaySpecies = {'CO2', 'CO', 'H2O', 'H2', 'O2', 'H','OH','O',...
                      'CH4','C2H4','CH3','HCO','CH','Cbgrb'};
    tolMoles = 1e-14;

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
    mixArray = setProperties(mix, 'temperature', 300, 'pressure', 1, 'equivalenceRatio', 0.5:0.01:4);
    
    % Initialize solver
    solver = EquilibriumSolver('problemType', 'HP', 'tolMoles', tolMoles, 'FLAG_RESULTS', false);
    
    % Solve problem
    solver.solveArray(mixArray);
    
    if FLAG_BENCHMARK
        stackTrace = dbstack; functionname = stackTrace.name;
        metadata = combustiontoolbox.utils.BenchmarkMetadata(solver, mixArray, functionname);
        return
    end

    % Load results CEA 
    resultsCEA = data_CEA(filename, displaySpecies);
    
    % Plot molar fractions
    plotComposition(mixArray(1), mixArray, 'equivalenceRatio', 'Xi', 'displaySpecies', displaySpecies, 'mintol', 1e-14, 'validation', resultsCEA);
    
    % Plot molar fractions
    fig1 = plotComposition(mixArray(1), mixArray, 'equivalenceRatio', 'Xi', 'displaySpecies', displaySpecies, 'mintol', 1e-14, 'validation', resultsCEA);

    % Plot properties
    fig2 = plotProperties(repmat({'equivalenceRatio'}, 1, 9), mixArray, {'T', 'rho', 'h', 'e', 'g', 'cp', 's', 'gamma_s', 'sound'}, mixArray, 'basis', {[], [], 'mi', 'mi', 'mi', 'mi', 'mi', [], []}, 'validation', resultsCEA);

    % Save plots
    if ~FLAG_EXPORT
        return
    end

    folderpath = fullfile(pwd, 'validations', 'figures');

    % Check if folder exists, otherwise create it
    if ~exist(folderpath, 'dir')
        mkdir(folderpath);
    end

    stackTrace = dbstack;
    filename = stackTrace.name;
    saveas(fig1, fullfile(folderpath, strcat(filename, '_molar')), 'svg');
    saveas(fig2, fullfile(folderpath, strcat(filename, '_properties')), 'svg');
end