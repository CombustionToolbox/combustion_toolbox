function run_validation_SHOCK_POLAR_SDToolbox_2
    % Run test validation_SHOCK_POLAR_SDToolbox_2
    % Contrasted with: Caltech's SD Toolbox and CANTERA
    % Problem type: Shock polar
    % Temperature [K]   = 300
    % Pressure    [bar] = 1.01325
    % Initial mixture: AIR (79% N2 + 21% O2)
    % List of species considered: list_species('AIR_IONS')

    % Import packages
    import combustiontoolbox.databases.NasaDatabase
    import combustiontoolbox.core.*
    import combustiontoolbox.shockdetonation.*
    import combustiontoolbox.utils.display.*

    % Benchmark?
    FLAG_BENCHMARK = false;

    % Definitions
    prefixDataName = 'shock_polar_equilwithions_ar_SDToolbox';

    for i = 6:-1:1
        filename{i} = sprintf('%s%d.out', prefixDataName, i);
    end

    listSpecies = {'O2', 'O', 'N2', 'N', 'NO', 'eminus', 'Nplus', 'NOplus',...
          'N2plus', 'Oplus', 'O2plus', 'Ar', 'Arplus'};
    M1 = [2, 3, 5, 14];
    numPointsPolar = 300;
    filename = 'shock_polar_equilwithions_SDToolbox';
    
    % Get Nasa database
    DB = NasaDatabase('FLAG_BENCHMARK', FLAG_BENCHMARK);
    
    % Define chemical system
    system = ChemicalSystem(DB, listSpecies);
    
    % Initialize mixture
    mix = Mixture(system);
    
    % Define chemical state
    set(mix, {'N2', 'O2', 'Ar'}, [78, 21, 1] / 21);
    
    % Define properties
    mixArray1 = setProperties(mix, 'temperature', 300, 'pressure', 1, 'M1', M1);
    
    % Initialize solver
    solver = ShockSolver('problemType', 'SHOCK_POLAR', 'numPointsPolar', numPointsPolar, 'FLAG_RESULTS', false);
    
    % Solve problem
    [mixArray1, mixArray2] = solver.solveArray(mixArray1);
    
    if FLAG_BENCHMARK
        return
    end

    % Load results SDToolbox 
    resultsSDToolbox = load_struct(filename, 'data');

    % Plot polars
    [ax1, ax2] = plot_validation_shock_polar_SDToolbox(mixArray1, mixArray2, resultsSDToolbox, PlotConfig());

    % Save plots
    folderpath = strcat(fullfile(pwd, 'validations', 'figures'));
    stack_trace = dbstack;
    filename = stack_trace.name;
    saveas(ax1, fullfile(folderpath, strcat(filename, '_properties_1')), 'svg');
    saveas(ax2, fullfile(folderpath, strcat(filename, '_properties_2')), 'svg');
end