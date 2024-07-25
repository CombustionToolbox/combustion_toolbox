function [max_rel_error_moles, max_rel_error_prop] = run_test_DET_CEA_2(value, database)
    % Run test_DET_CEA_1:
    % Contrasted with: NASA's Chemical Equilibrium with Applications software
    % Problem type: Chapman-Jouguet Detonation
    % Temperature [K]   = 300;
    % Pressure    [bar] = 1;
    % Equivalence ratio [-] = 0.5:0.01:4
    % Initial mixture: Fuel + O2
    % List of species considered: All possible products
    
    % Import packages
    import combustiontoolbox.databases.NasaDatabase
    import combustiontoolbox.core.*
    import combustiontoolbox.shockdetonation.DetonationSolver
    import combustiontoolbox.utils.display.*
    
    % If equivalence ratio == 3.9, return directly. CEA does not converge
    % for that value.
    if value == 3.9
        max_rel_error_moles = 0;
        max_rel_error_prop = 0;
        return
    end

    % Definitions
    fuel = 'C2H2_acetylene';
    prefixDataName = 'C2H2_O2_detonations';

    for i = 3:-1:1
        filename{i} = sprintf('%s%d.out', prefixDataName, i);
    end

    displaySpecies = {'CO2', 'CO', 'H2O', 'H2', 'O2', 'H','OH','O',...
                      'CH4','C2H4','CH3','HCO','CH','Cbgrb'};
    tolMoles = 1e-18;
    propertiesY = {'T', 'p', 'h', 'e', 'g', 's', 'cp', 'cv', 'gamma_s', 'dVdT_p', 'dVdp_T', 'sound', 'W', 'uShock'};
    propertiesX = repmat({'equivalenceRatio'}, 1, length(propertiesY));
    
    % Define chemical system
    system = ChemicalSystem(database);
    
    % Initialize mixture
    mix = Mixture(system);
    
    % Define chemical state
    set(mix, {fuel}, 'fuel', 1);
    set(mix, {'O2'}, 'oxidizer', 1);
    
    % Define properties
    mixArray1 = setProperties(mix, 'temperature', 300, 'pressure', 1, 'equivalenceRatio',  value);
    
    % Initialize solver
    solver = DetonationSolver('problemType', 'DET', 'FLAG_RESULTS', true, 'tolMoles', tolMoles);
    
    % Solve problem
    [mixArray1, mixArray2] = solver.solveArray(mixArray1);
    
    % Prepare data
    for i = 1:length(mixArray2)
        mixArray2(i).uShock = mixArray1(i).u;
    end

    % Load results CEA 
    resultsCEA = data_CEA(filename, displaySpecies);
    resultsCEA.uShock = resultsCEA.u_preshock;

    % Compute error: molar fractions
    max_rel_error_moles = compute_error_moles_CEA(mixArray2, resultsCEA, 'equivalenceRatio', value, 'Xi', displaySpecies, tolMoles);
    
    % Compute error: properties mixture 2
    max_rel_error_prop = compute_error_prop_CEA(mixArray2, resultsCEA, propertiesX, value, propertiesY, 'mix2');
end