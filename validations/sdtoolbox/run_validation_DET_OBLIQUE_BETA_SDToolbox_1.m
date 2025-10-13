function metadata = run_validation_DET_OBLIQUE_BETA_SDToolbox_1(varargin)
    % Run test validation_DET_OBLIQUE_BETA_SDToolbox_1
    % Contrasted with: Caltech's SD Toolbox and CANTERA
    % Problem type: Oblique detonation wave given wave angle
    % Temperature [K]   = 300
    % Pressure    [bar] = 1.01325
    % Initial mixture: CH4 + AIR (79% N2 + 21% O2)
    % Overdrive factor = 4
    % Beta angle = 15:1:90
    % List of species considered: 
    
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
    fuel = 'CH4';
    prefixDataName = 'airNASA_incident_shocks_ionization';

    for i = 6:-1:1
        filename{i} = sprintf('%s%d.out', prefixDataName, i);
    end

    listSpecies = 'Soot Formation';
    equivalenceRatio = 1;
    driveFactor = 4;
    beta = 15:1:90;
    filename = 'shock_polar_equilwithions_SDToolbox';
    
    % Get Nasa database
    DB = NasaDatabase('FLAG_BENCHMARK', FLAG_BENCHMARK);
    
    % Define chemical system
    system = ChemicalSystem(DB, listSpecies);
    
    % Initialize mixture
    mix = Mixture(system);
    
    % Define chemical state
    set(mix, {fuel}, 'fuel', 1);
    set(mix, {'N2', 'O2'}, 'oxidizer', [79, 21] / 21);
    
    % Define properties
    mixArray1 = setProperties(mix, 'temperature', 300, 'pressure', 1, 'equivalenceRatio', equivalenceRatio, 'driveFactor', driveFactor, 'beta', beta);
    
    % Initialize solver
    solver = DetonationSolver('problemType', 'DET_OBLIQUE', 'FLAG_RESULTS', false);
    
    % Solve problem
    [mixArray1, mixArray2_1, mixArray2_2] = solver.solveArray(mixArray1);
    
    if FLAG_BENCHMARK
        stackTrace = dbstack; functionname = stackTrace.name;
        metadata = combustiontoolbox.utils.BenchmarkMetadata(solver, mixArray1, functionname);
        return
    end

    % % Load results SDToolbox 
    % resultsSDToolbox = load_struct(filename, 'data');

    % % Plot polars
    % [ax1, ax2] = plot_validation_shock_polar_SDToolbox(mixArray1, mixArray2, resultsSDToolbox, PlotConfig());

    % % Save plots
    % if ~FLAG_EXPORT
    %     return
    % end

    % folderpath = fullfile(pwd, 'validations', 'figures');

    % % Check if folder exists, otherwise create it
    % if ~exist(folderpath, 'dir')
    %     mkdir(folderpath);
    % end

    % stackTrace = dbstack;
    % filename = stackTrace.name;
    % saveas(ax1, fullfile(folderpath, strcat(filename, '_properties_1')), 'svg');
    % saveas(ax2, fullfile(folderpath, strcat(filename, '_properties_2')), 'svg');
end