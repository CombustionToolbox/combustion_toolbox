function run_validation_DET_CEA_2
    % Run test validation_DET_CEA_2:
    % Contrasted with: NASA's Chemical Equilibrium with Applications software
    % Problem type: Chapman-Jouguet Detonation
    % Temperature [K]   = 300;
    % Pressure    [bar] = 1;
    % Equivalence ratio [-] = 0.5:0.01:4
    % Initial mixture: C2H2_acetylene + O2
    % List of species considered: All possible products
    
    % Import packages
    import combustiontoolbox.databases.NasaDatabase
    import combustiontoolbox.core.*
    import combustiontoolbox.shockdetonation.*
    import combustiontoolbox.utils.display.*

    % Benchmark?
    FLAG_BENCHMARK = false;

    % Definitions
    fuel = 'C2H2_acetylene';
    prefixDataName = 'C2H2_O2_detonations';

    for i = 3:-1:1
        filename{i} = sprintf('%s%d.out', prefixDataName, i);
    end

    displaySpecies = {'CO2', 'CO', 'H2O', 'H2', 'O2', 'H','OH','O',...
                      'CH4','C2H4','CH3','HCO','CH','Cbgrb'};
    
    % Get Nasa database
    DB = NasaDatabase('FLAG_BENCHMARK', FLAG_BENCHMARK);
    
    % Define chemical system
    system = ChemicalSystem(DB);
    
    % Initialize mixture
    mix = Mixture(system);
    
    % Define chemical state
    set(mix, {fuel}, 'fuel', 1);
    set(mix, {'O2'}, 'oxidizer', 1);
    
    % Define properties
    mixArray1 = setProperties(mix, 'temperature', 300, 'pressure', 1, 'equivalenceRatio',  0.5:0.01:4);
    
    % Initialize solver
    solver = DetonationSolver('problemType', 'DET', 'FLAG_RESULTS', false);
    
    % Solve problem
    [mixArray1, mixArray2] = solver.solveArray(mixArray1);
    
    if FLAG_BENCHMARK
        return
    end
    
    % Prepare data
    for i = 1:length(mixArray2)
        mixArray2(i).W = mixArray2(i).W * 1e3; % [g/mol]
        mixArray2(i).uShock = mixArray1(i).u;
    end

    % Load results CEA 
    resultsCEA = data_CEA(filename, displaySpecies);
    resultsCEA.uShock = resultsCEA.u_preshock;
    
    % Plot molar fractions
    fig1 = plotComposition(mixArray2(1), mixArray1, 'equivalenceRatio', 'Xi', 'mintol', 1e-14, 'y_var', mixArray2, 'validation', resultsCEA, 'display_species', displaySpecies);

    % Properties mixture 2 - 1
    fig2 = plotProperties(repmat({'equivalenceRatio'}, 1, 8), [mixArray1.equivalenceRatio], {'T', 'p', 'rho', 'h', 'e', 'g', 's', 'gamma_s'}, mixArray2, 'basis', {[], [], [], 'mi', 'mi', 'mi', 'mi', []}, 'validation', resultsCEA);

    % Properties mixture 2 - 2
    fig3 = plotProperties(repmat({'equivalenceRatio'}, 1, 7), [mixArray1.equivalenceRatio], {'cp', 'cv', 'dVdT_p', 'dVdp_T', 'sound', 'W', 'uShock'}, mixArray2, 'basis', {'mi', 'mi', [], [], [], [], []}, 'validation', resultsCEA);

    % Save plots
    folderpath = strcat(pwd,'\validations\figures\');
    stack_trace = dbstack;
    filename = stack_trace.name;
    saveas(fig1, strcat(folderpath, filename, '_molar'), 'svg');
    saveas(fig2, strcat(folderpath, filename, '_properties_1'), 'svg');
    saveas(fig3, strcat(folderpath, filename, '_properties_2'), 'svg');
end