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
    plotComposition(mixArray(1), mixArray, 'equivalenceRatio', 'Xi', 'displaySpecies', displaySpecies, 'mintol', 1e-14, 'validation', resultsCEA);
end

% % Display validation (plot)
% % * Molar fractions
% [~, fig1] = plot_molar_fractions(results_CT, results_CT.PS.strP, 'phi', 'Xi', 'validation', results_CEA, 'display_species', display_species);
% % * Properties mixture 2
% fig2 = plot_properties_validation(results_CT, results_CEA, {'phi', 'phi', 'phi', 'phi', 'phi', 'phi', 'phi', 'phi'}, {'rho', 'h', 'e', 'g', 'S', 'cP', 'cV', 'gamma_s'}, 'mix2');
% % Save plots
% folderpath = strcat(pwd,'\Validations\Figures\');
% stack_trace = dbstack;
% filename = stack_trace.name;
% saveas(fig1, strcat(folderpath, filename, '_molar'), 'svg');
% saveas(fig2, strcat(folderpath, filename, '_properties'), 'svg');