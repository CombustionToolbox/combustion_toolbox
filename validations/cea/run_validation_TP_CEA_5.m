function run_validation_TP_CEA_5
    % Run test validation_TP_CEA_5:
    % Contrasted with: NASA's Chemical Equilibrium with Applications software
    % Problem type: Equilibrium composition at defined T and p
    % Temperature [K]   = 2000:100:20000;
    % Pressure    [bar] = 1;
    % Equivalence ratio [-] = 1.1
    % Initial mixture: C6H6 + AIR_IDEAL (79% N2 + 21% O2)
    % List of species considered: list_species('NASA ALL IONS')

    % Import packages
    import combustiontoolbox.databases.NasaDatabase
    import combustiontoolbox.core.*
    import combustiontoolbox.equilibrium.*
    import combustiontoolbox.utils.display.*
    
    % Benchmark?
    FLAG_BENCHMARK = false;

    % Definitions
    fuel = 'C6H6';
    prefixDataName = fuel;
    filename = {strcat(prefixDataName, '_air_T_TP1.out'), strcat(prefixDataName, '_air_T_TP2.out')};
    listSpecies = 'NASA ALL IONS';
    displaySpecies = {'CO2', 'CO', 'H2O', 'H2', 'O2', 'N2', 'He', 'Ar',...
                      'HCN','H','OH','O','CN','NH3','CH4','C2H4','CH3',...
                      'NO','HCO','NH2','NH','N','CH','Cbgrb'};
    tolMoles = 1e-18;

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
    mixArray = setProperties(mix, 'temperature', 2000:100:20000, 'pressure', 1, 'equivalenceRatio', 1.1);
    
    % Initialize solver
    solver = EquilibriumSolver('problemType', 'TP', 'tolMoles', tolMoles, 'FLAG_RESULTS', false);
    
    % Solve problem
    solver.solveArray(mixArray);
    
    if FLAG_BENCHMARK
        return
    end

    % Load results CEA 
    resultsCEA = data_CEA(filename, displaySpecies);
    
    % Plot molar fractions
    fig1 = plotComposition(mixArray(1), mixArray, 'equivalenceRatio', 'Xi', 'displaySpecies', displaySpecies, 'mintol', 1e-3, 'validation', resultsCEA);

    % Properties mixture
    fig2 = plotProperties(repmat({'equivalenceRatio'}, 1, 10), mixArray, {'rho', 'h', 'e', 'g', 's', 'cp', 'cv', 'gamma_s', 'dVdp_T', 'dVdT_p'}, mixArray, 'basis', {[], 'mi', 'mi', 'mi', 'mi', 'mi', 'mi', [], [], []}, 'validation', resultsCEA);

    % Save plots
    folderpath = fullfile(pwd, 'validations', 'figures');
    stack_trace = dbstack;
    filename = stack_trace.name;
    saveas(fig1, fullfile(folderpath, strcat(filename, '_molar')), 'svg');
    saveas(fig2, fullfile(folderpath, strcat(filename, '_properties')), 'svg');
end