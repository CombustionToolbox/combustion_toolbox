 function metadata = run_validation_SHOCK_IONIZATION_CEA_2(varargin)
    % Run test validation_SHOCK_IONIZATION_CEA_2:
    % Contrasted with: NASA's Chemical Equilibrium with Applications software
    % Problem type: Planar incident shock wave
    % Temperature [K]   = 300
    % Pressure    [bar] = 1
    % Incident velocity [m/s] = [~308, 13000]
    % Initial mixture: AIR (78.084% N2 + 20.9476% O2 + 0.9365% Ar + 0.0319% CO2)
    % List of species considered: All + ions (see method findProducts from ChemicalSystem class)
    
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
    prefixDataName = 'airNASA_incident_shocks_ionization';

    for i = 6:-1:1
        filename{i} = sprintf('%s%d.out', prefixDataName, i);
    end

    listSpecies = [];
    displaySpecies = {'eminus','Ar','Arplus','C','Cplus','Cminus','CN','CNplus','CNminus',...
                      'CNN','CO','COplus','CO2','CO2plus',...
                      'N','Nplus','Nminus','NCO','NO','NOplus','NO2','NO2minus',...
                      'N2','N2plus','N2minus','NCN','N2O','N2Oplus','N3',...
                      'O','Oplus','Ominus','O2','O2plus','O2minus','O3'};
    u1 = logspace(2, 5, 500); u1 = u1(u1<13000); u1 = u1(u1>=357);
    
    % Get Nasa database
    DB = NasaDatabase('FLAG_BENCHMARK', FLAG_BENCHMARK);
    
    % Define chemical system
    system = ChemicalSystem(DB);
    system.FLAG_ION = true;
    
    % Initialize mixture
    mix = Mixture(system);
    
    % Define chemical state
    set(mix, {'N2', 'O2', 'Ar', 'CO2'}, [78.084, 20.9476, 0.9365, 0.0319] / 20.9476);
    
    % Define properties
    mixArray1 = setProperties(mix, 'temperature', 300, 'pressure', 1, 'u1', u1);
    
    % Initialize solver
    solver = ShockSolver('problemType', 'SHOCK_I', 'FLAG_RESULTS', false);
    
    % Solve problem
    [mixArray1, mixArray2] = solver.solveArray(mixArray1);
    
    if FLAG_BENCHMARK
        stackTrace = dbstack; functionname = stackTrace.name;
        metadata = combustiontoolbox.utils.BenchmarkMetadata(solver, mixArray1, functionname);
        return
    end
    
    % Prepare data
    for i = 1:length(mixArray2)
        mixArray2(i).W = mixArray2(i).W * 1e3; % [g/mol]
    end

    % Load results CEA 
    resultsCEA = data_CEA(filename, displaySpecies);
    
    % Plot molar fractions
    fig1 = plotComposition(mixArray2(1), mixArray1, 'u', 'Xi', 'mintol', 1e-5, 'y_var', mixArray2);

    % Properties mixture 1
    fig2 = plotProperties(repmat({'u_preshock'}, 1, 8), [mixArray1.u], {'T', 'p', 'rho', 'h', 'e', 'g', 's', 'gamma_s'}, mixArray1, 'basis', {[], [], [], 'mi', 'mi', 'mi', 'mi', []}, 'validation', resultsCEA.mix1);

    % Properties mixture 2 - 1
    fig3 = plotProperties(repmat({'u_preshock'}, 1, 8), [mixArray1.u], {'T', 'p', 'rho', 'h', 'e', 'g', 's', 'gamma_s'}, mixArray2, 'basis', {[], [], [], 'mi', 'mi', 'mi', 'mi', []}, 'validation', resultsCEA.mix2);

    % Properties mixture 2 - 2
    fig4 = plotProperties(repmat({'u_preshock'}, 1, 6), [mixArray1.u], {'cp', 'cv', 'dVdT_p', 'dVdp_T', 'sound', 'W'}, mixArray2, 'basis', {'mi', 'mi', [], [], [], []}, 'validation', resultsCEA.mix2);

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
    saveas(fig2, fullfile(folderpath, strcat(filename, '_properties_1')), 'svg');
    saveas(fig3, fullfile(folderpath, strcat(filename, '_properties_2')), 'svg');
    saveas(fig4, fullfile(folderpath, strcat(filename, '_properties_3')), 'svg');
 end