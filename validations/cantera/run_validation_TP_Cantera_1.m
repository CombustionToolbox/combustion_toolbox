function metadata = run_validation_TP_Cantera_1(varargin)
    % Run test validation_TP_Cantera_1:
    % Contrasted with: Cantera
    % Problem type: Frozen composition at defined T and p
    % Temperature [K]   = 200;
    % Pressure    [bar] = logspace(0, 3, 100);
    % Initial mixture: AIR_IDEAL (79% N2 + 21% O2)
    % Equation of state: Peng-Robinson
    
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
    prefixDataName = 'air';
    filename = strcat('Cantera_', prefixDataName, '_TP1.txt');
    listSpecies = [];

    % Get Nasa database
    DB = NasaDatabase('FLAG_BENCHMARK', FLAG_BENCHMARK);
    
    % Add Peng-Robinson parameters
    DB.species = addPengRobinsonProperties(DB.species);

    % Define chemical system
    system = ChemicalSystem(DB, listSpecies);
    
    % Define equation of state
    eos = EquationStatePengRobinson();

    % Initialize mixture
    mix = Mixture(system, 'eos', eos);
    
    % Define chemical state
    set(mix, {'N2', 'O2'}, [79, 21] / 21);
    
    % Define properties
    mixArray = setProperties(mix, 'temperature', 200, 'pressure', logspace(0, 2, 100));

    % Define caloric gas model
    caloricGasModel = CaloricGasModel.thermallyPerfect;

    % Initialize solver
    solver = EquilibriumSolver('problemType', 'TP', 'caloricGasModel', caloricGasModel, 'FLAG_RESULTS', false);

    % Solve problem
    solver.solveArray(mixArray);

    if FLAG_BENCHMARK
        stackTrace = dbstack; functionname = stackTrace.name;
        metadata = combustiontoolbox.utils.BenchmarkMetadata(solver, mixArray, functionname);
        return
    end

    % Load results Cantera 
    resultsCantera = readtable(filename, 'PreserveVariableNames', true);

    % Properties mixture
    properties = {'rho', 'h', 'e', 'g', 's', 'cp', 'cv', 'gamma_s', 'dVdp_T', 'dVdT_p', 'dPdV_T', 'dPdT_V'};
    basis = {[], 'mi', 'mi', 'mi', 'mi', 'mi', 'mi', [], [], [], [], []};

    fig1 = plotProperties(repmat({'p'}, 1, length(properties)), mixArray, properties, mixArray, 'basis', basis, 'validation', resultsCantera);

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
    saveas(fig1, fullfile(folderpath, strcat(filename, '_properties')), 'svg');
end