function run_validation_TP_CEA_6(varargin)
    % Run test validation_TP_CEA_6:
    % Contrasted with: NASA's Chemical Equilibrium with Applications software
    % Problem type: Equilibrium composition at defined T and p
    % Temperature [K]   = 2000:50:5000;
    % Pressure    [bar] = 1.01325;
    % Initial mixture: Si + 9 C6H5OH_phenol
    % List of species considered: All (see method findProducts from ChemicalSystem class)
    %
    % This example is obtained from [1]
    % 
    % [1] Scoggins, J. B., & Magin, T. E. (2015). Gibbs function continuation
    %     for linearly constrained multiphase equilibria. Combustion and Flame,
    %     162(12), 4514-4522.
  
    % Import packages
    import combustiontoolbox.databases.NasaDatabase
    import combustiontoolbox.core.*
    import combustiontoolbox.equilibrium.*
    import combustiontoolbox.utils.display.*
    
    % Benchmark?
    FLAG_BENCHMARK = false;

    % Definitions
    displaySpecies = {'H2', 'H', 'CH4', 'C2H2_acetylene', 'SiC2',...
        'Si', 'C', 'C2H', 'C3', 'CO', 'CO2', 'Cbgrb', 'SiCbbb', 'H2O',...
        'SiO2ba_qzb','SiO2bb_qzb','SiO2bb_crtb','H2ObLb','H2Obcrb'};

    % Get Nasa database
    DB = NasaDatabase('FLAG_BENCHMARK', FLAG_BENCHMARK);
    
    % Define chemical system
    system = ChemicalSystem(DB);
    
    % Initialize mixture
    mix = Mixture(system);
    
    % Define chemical state
    set(mix, {'Si', 'C6H5OH_phenol'}, [1, 9]);
    
    % Define properties
    mixArray = setProperties(mix, 'temperature', 200:10:5000, 'pressure', 1 * 1.01325);
    
    % Initialize solver
    solver = EquilibriumSolver('problemType', 'TP', 'FLAG_RESULTS', false);
    
    % Solve problem
    solver.solveArray(mixArray);
    
    if FLAG_BENCHMARK
        return
    end

    % Load results CEA 
    prefixDataName = 'C6H5OH_phenol_and_Si';
    filename = {strcat(prefixDataName, '_TP1.out')};
    resultsCEA = data_CEA(filename, displaySpecies);

    % Plot molar fractions
    fig1 = plotComposition(mixArray(1), mixArray, 'T', 'Xi', 'displaySpecies', displaySpecies, 'mintol', 1e-3, 'validation', resultsCEA);

    % Properties mixture
    fig2 = plotProperties(repmat({'T'}, 1, 10), mixArray, {'rho', 'h', 'e', 'g', 's', 'cp', 'cv', 'gamma_s', 'dVdp_T', 'dVdT_p'}, mixArray, 'basis', {[], 'mi', 'mi', 'mi', 'mi', 'mi', 'mi', [], [], []}, 'validation', resultsCEA);

    % Save plots
    folderpath = fullfile(pwd, 'validations', 'figures');
    stack_trace = dbstack;
    filename = stack_trace.name;
    saveas(fig1, fullfile(folderpath, strcat(filename, '_molar')), 'svg');
    saveas(fig2, fullfile(folderpath, strcat(filename, '_properties')), 'svg');
end