function run_validation_TP_TEA_3
    % VALIDATION: TP_TEA_3
    %
    % Compute equilibrium composition at defined temperature and pressure.
    % Reproduce Fig.6 (b) "Thermochemical equilibrium vertical distributions
    % with a metallicity 10 of WASP-43b assuming the T-P profile in Fig 4"
    % [TEA by Jasmina Blecic].
    %
    % URL RESULTS TEA:
    % https://github.com/dzesmin/RRC-BlecicEtal-2015a-ApJS-TEA/tree/master/Fig6/WASP43b-10xsolar
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
    load Validation_TP_TEA_3 Pressure Temp results_TEA

    % Definitions
    listSpecies = {'C2H2_acetylene', 'C2H4', 'C', 'CH4', 'CO2', 'CO', 'H2',...
          'H2O', 'H2S', 'H', 'HCN', 'He', 'HS_M', 'N2', 'N', 'NH3',...
          'O', 'S', 'HSO_M', 'HSO2_M', 'HSO3_M', 'HS2_M', 'S2',...
          'S2O_M', 'S_OH_M', 'OH'};

    displaySpecies = {'C2H2_acetylene', 'C2H4', 'C', 'CH4', 'CO2', 'CO', 'H2',...
      'H2O', 'H2S', 'H', 'HCN', 'He', 'HS_M', 'N2', 'N', 'NH3',...
      'O', 'S'};

    species = {'H', 'He', 'C', 'N', 'O', 'S'};
    metallicity = 10;
    
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
    mixArray = setProperties(mix, 'temperature', Temp, 'pressure', Pressure);
    
    % Initialize solver
    solver = EquilibriumSolver('problemType', 'TP', 'FLAG_RESULTS', false, 'tolMoles', 1e-32);
    
    % Solve problem
    solver.solveArray(mixArray);
    
    if FLAG_BENCHMARK
        return
    end

    % Plot molar fractions
    fig1 = plotComposition(mixArray(1), mixArray, 'Xi', 'p', 'displaySpecies', displaySpecies, 'mintol', 1e-21, 'nfrec', 3, 'ydir', 'reverse', 'xscale', 'log', 'validation', results_TEA);

    % Save plots
    folderpath = fullfile(pwd, 'validations', 'figures');
    stack_trace = dbstack;
    filename = stack_trace.name;
    saveas(fig1, fullfile(folderpath, strcat(filename, '_molar')), 'svg');
end
