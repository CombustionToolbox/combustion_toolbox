function metadata = run_validation_SHOCK_PRANDTL_MEYER_SDToolbox_2(varargin)
    % Run test validation_SHOCK_PRANDTL_MEYER_SDToolbox_2
    % Contrasted with: Caltech's SD Toolbox and CANTERA
    % Problem type: Shock Prandtl-Meyer expansion (frozen)
    % Temperature [K]   = 3000
    % Pressure    [bar] = 1
    % Initial mixture: AIR (79% N2 + 21% O2)
    % Free-stream Mach number = 1
    % Wave angle = 0:80
    % List of species considered: list_species('AIR_IONS')
    
    % Import packages
    import combustiontoolbox.databases.NasaDatabase
    import combustiontoolbox.core.*
    import combustiontoolbox.shockdetonation.*
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
    filename = 'shock_prandtlmeyer_frozenwithions_SDToolbox1';
    listSpecies = {'O2', 'O', 'N2', 'N', 'NO', 'eminus', 'Nplus', 'NOplus',...
          'N2plus', 'Oplus', 'O2plus'};
    M1 = 1;
    theta = 0:1:80;
    numPointsPrandtlMeyer = 100;
    caloricGasModel = CaloricGasModel.thermallyPerfect;
    
    % Get Nasa database
    DB = NasaDatabase('FLAG_BENCHMARK', FLAG_BENCHMARK);
    
    % Define chemical system
    system = ChemicalSystem(DB, listSpecies);
    
    % Initialize mixture
    mix = Mixture(system);
    
    % Define chemical state
    set(mix, {'N2', 'O2'}, [79, 21] / 21);
    
    % Define properties
    mixArray1 = setProperties(mix, 'temperature', 3000, 'pressure', 1, 'M1', M1, 'theta', theta);
    
    % Initialize solver
    solver = ShockSolver('problemType', 'SHOCK_PRANDTL_MEYER', 'numPointsPrandtlMeyer', numPointsPrandtlMeyer, 'caloricGasModel', caloricGasModel, 'FLAG_RESULTS', false);
    
    % Solve problem
    [mixArray1, mixArray2] = solver.solveArray(mixArray1);
    
    if FLAG_BENCHMARK
        stackTrace = dbstack; functionname = stackTrace.name;
        metadata = combustiontoolbox.utils.BenchmarkMetadata(solver, mixArray1, functionname);
        return
    end

    % Load results SDToolbox 
    resultsSDToolbox = load_struct(filename, 'data');

    % Plot expansion path
    if isscalar(theta)
        mixArray2(end).polar(1).rangeName = 'theta';
        ax1 = plotProperties(repmat({'theta'}, 1, 5), [mixArray2(end).polar.theta], {'T', 'p', 'rho', 'mach', 'gamma_f'}, mixArray2(end).polar, 'basis', {[], [], [], [], []}, 'validation', resultsSDToolbox);
    else
        ax1 = plotProperties(repmat({'theta'}, 1, 5), [mixArray2.theta], {'T', 'p', 'rho', 'mach', 'gamma_f'}, mixArray2, 'basis', {[], [], [], [], []}, 'validation', resultsSDToolbox);
    end

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
    saveas(ax1, fullfile(folderpath, strcat(filename, '_properties_1')), 'svg');
end