function [max_rel_error_moles, max_rel_error_prop] = run_test_HP_CEA_1(value, database)
    % Run test_HP_CEA_1:
    % Contrasted with: NASA's Chemical Equilibrium with Applications software
    % Problem type: Adiabatic T and composition at constant p
    % Temperature [K]   = 300;
    % Pressure    [bar] = 1;
    % Equivalence ratio [-] = value
    % Initial mixture: Fuel + AIR_IDEAL (79% N2 + 21% O2)
    % List of species considered: list_species('Soot Formation Extended')
    
    % Import packages
    import combustiontoolbox.databases.NasaDatabase
    import combustiontoolbox.core.*
    import combustiontoolbox.equilibrium.*
    import combustiontoolbox.utils.display.*

    % Definitions
    fuel = 'C2H2_acetylene';
    prefixDataName = 'C2H2';
    filename = {strcat(prefixDataName, '_air1_HP.out'), strcat(prefixDataName, '_air1_HP2.out'), strcat(prefixDataName, '_air1_HP3.out')};
    displaySpecies = {'CO2', 'CO', 'H2O', 'H2', 'O2', 'N2',...
                      'HCN', 'H', 'OH', 'O', 'CN', 'NH3', 'CH4', 'C2H4', ...
                      'CH3', 'NO', 'HCO', 'NH2', 'NH', 'N', 'CH', 'Cbgrb'};
    tolMoles = 1e-18;
    propertiesY = {'T', 'p', 'h', 'e', 'g', 's', 'cp', 'cv', 'gamma_s', 'dVdT_p', 'dVdp_T', 'sound', 'W'};
    propertiesX = repmat({'equivalenceRatio'}, 1, length(propertiesY));

    % Define chemical system
    system = ChemicalSystem(database);
    
    % Initialize mixture
    mix = Mixture(system);
    
    % Define chemical state
    set(mix, {fuel}, 'fuel', 1);
    set(mix, {'N2', 'O2'}, 'oxidizer', [79, 21] / 21);
    
    % Define properties
    mixArray = setProperties(mix, 'temperature', 300, 'pressure', 1, 'equivalenceRatio', value);
    
    % Initialize solver
    solver = EquilibriumSolver('problemType', 'HP', 'FLAG_RESULTS', false, 'tolMoles', tolMoles);
    
    % Solve problem
    solver.solveArray(mixArray);
    
    % Load results CEA 
    resultsCEA = data_CEA(filename, displaySpecies);

    % Compute error: molar fractions
    max_rel_error_moles = compute_error_moles_CEA(mixArray, resultsCEA, 'equivalenceRatio', value, 'Xi', displaySpecies, tolMoles);
    
    % Compute error: properties mixture 2
    max_rel_error_prop = compute_error_prop_CEA(mixArray, resultsCEA, propertiesX, value, propertiesY, 'mix2');
end