function run_validation_TP_TEA_1
    % VALIDATION: TP_TEA_1
    %
    % Compute equilibrium composition at defined temperature and pressure.
    % Reproduce the example case of TEA by Jasmina Blecic.
    % URL RESULTS TEA:
    % https://github.com/dzesmin/TEA/tree/master/doc/examples/quick_example/results 
    %   
    %
    % @author: Alberto Cuadra Lara
    %          Postdoctoral researcher - Group Fluid Mechanics
    %          Universidad Carlos III de Madrid
    %                 
    % Last update Jun 06 2024
    
    % Import packages
    import combustiontoolbox.databases.SolarAbundances
    import combustiontoolbox.databases.NasaDatabase
    import combustiontoolbox.core.*
    import combustiontoolbox.equilibrium.*
    import combustiontoolbox.utils.display.*
    
    % Benchmark?
    FLAG_BENCHMARK = false;

    % Inputs
    load Validation_TP_TEA_1 Temp results_TEA

    % Definitions
    listSpecies = {'C', 'CH4', 'CO2', 'CO', 'H2', 'H', 'H2O', 'He', 'N2', 'N',...
          'NH3', 'O'};

    displaySpecies = listSpecies;

    species = {'H', 'He', 'C', 'N', 'O'};
    metallicity = 1;
    
    % Get initial composition from solar abundances
    DB_solar = SolarAbundances();
    moles = DB_solar.abundances2moles(species, metallicity);
    
    % Get Nasa database
    DB = NasaDatabase('FLAG_BENCHMARK', true);
    
    % Define chemical system
    system = ChemicalSystem(DB, listSpecies);
    
    % Initialize mixture
    mix = Mixture(system);
    
    % Define chemical state
    set(mix, species, moles);
    
    % Define properties
    mixArray = setProperties(mix, 'temperature', linspace(Temp(1), Temp(end), 300), 'pressure', logspace(-5, 2, 300));
    
    % Initialize solver
    solver = EquilibriumSolver('problemType', 'TP', 'FLAG_RESULTS', false, 'tolMoles', 1e-32);
    
    % Solve problem
    solver.solveArray(mixArray);
    
    if FLAG_BENCHMARK
        return
    end

    % Plot molar fractions
    fig1 = plotComposition(mixArray(1), mixArray, 'T', 'Xi', 'displaySpecies', displaySpecies, 'mintol', 1e-14, 'nfrec', 3, 'validation', results_TEA);

    % Save plots
    folderpath = fullfile(pwd, 'validations', 'figures');
    stack_trace = dbstack;
    filename = stack_trace.name;
    saveas(fig1, fullfile(folderpath, strcat(filename, '_molar')), 'svg');
end
