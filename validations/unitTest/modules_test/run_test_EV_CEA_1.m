function [max_rel_error_moles, max_rel_error_prop] = run_test_EV_CEA_1(value, database)
    % Run test_EV_CEA_1:
    % Contrasted with: NASA's Chemical Equilibrium with Applications software
    % Problem type: Adiabatic T and composition at constant volume
    % Temperature [K]   = 300;
    % Specific volume [m3/kg] = 1;
    % Equivalence ratio [-] = value
    % Initial mixture: Fuel + AIR (78.084% N2, 20.9476% O2, 0.9365% Ar, 0.0319% CO2)
    % List of species considered: list_species('Soot Formation Extended')
    
        % Import packages
    import combustiontoolbox.databases.NasaDatabase
    import combustiontoolbox.core.*
    import combustiontoolbox.equilibrium.*
    import combustiontoolbox.utils.display.*

    % Definitions
    fuel = 'CH4';
    prefixDataName = fuel;

    for i = 2:-1:1
        filename{i} = strcat(prefixDataName, '_air_EV', sprintf('%d', i), '.out');
    end

    displaySpecies = {'CO2', 'CO', 'H2O', 'H2', 'O2', 'N2',...
                      'HCN', 'H', 'OH', 'O', 'CN', 'NH3', 'CH4', 'C2H4', ...
                      'CH3', 'NO', 'HCO', 'NH2', 'NH', 'N', 'CH', 'Cbgrb'};
    tolMoles = 1e-18;
    propertiesY = {'T', 'p', 'h', 'e', 'g', 's', 'cp', 'cv', 'gamma_s', 'dVdT_p', 'dVdp_T', 'sound', 'W'};
    propertiesX = repmat({'equivalenceRatio'}, size(propertiesY));

    % Define chemical system
    system = ChemicalSystem(database);
    
    % Initialize mixture
    mix = Mixture(system);
    
    % Define chemical state
    set(mix, {fuel}, 'fuel', 1);
    set(mix, {'N2', 'O2', 'Ar', 'CO2'}, 'oxidizer', [78.084, 20.9476, 0.9365, 0.0319] / 20.9476);
    
    % Define properties
    mixArray = setProperties(mix, 'temperature', 300, 'volume', 1, 'equivalenceRatio', value);
    
    % Initialize solver
    solver = EquilibriumSolver('problemType', 'EV', 'FLAG_RESULTS', false, 'tolMoles', tolMoles);
    
    % Solve problem
    solver.solveArray(mixArray);
    
    % Load results CEA 
    results_CEA = data_CEA(filename, displaySpecies);

    % Compute error: molar fractions
    max_rel_error_moles = compute_error_moles_CEA(mixArray, results_CEA, 'equivalenceRatio', value, 'Xi', displaySpecies, tolMoles);
    
    % Compute error: properties mixture 2
    max_rel_error_prop = compute_error_prop_CEA(mixArray, results_CEA, propertiesX, value, propertiesY, 'mix2');
end