function [max_rel_error_moles, max_rel_error_prop] = run_test_DET_CEA_1(value, database)
    % Run test_DET_CEA_1:
    % Contrasted with: NASA's Chemical Equilibrium with Applications software
    % Problem type: Chapman-Jouguet Detonation
    % Temperature [K]   = 300;
    % Pressure    [bar] = 1;
    % Equivalence ratio [-] = 0.5:0.01:4
    % Initial mixture: Fuel + AIR_IDEAL (79% N2 + 21% O2)
    % List of species considered: list_species('Soot Formation Extended')
    
    % Import packages
    import combustiontoolbox.databases.NasaDatabase
    import combustiontoolbox.core.*
    import combustiontoolbox.shockdetonation.DetonationSolver
    import combustiontoolbox.utils.display.*
    
    % Definitions
    fuel = 'C2H2_acetylene';
    prefixDataName = 'C2H2_air1_detonations';

    for i = 3:-1:1
        filename{i} = sprintf('%s%d.out', prefixDataName, i);
    end

    listSpecies = 'Soot Formation Extended';
    displaySpecies = {'CO2', 'CO', 'H2O', 'H2', 'O2', 'N2',...
                      'HCN','H','OH','O','CN','NH3','CH4','C2H4','CH3',...
                      'NO','HCO','NH2','NH','N','CH','Cbgrb'};
    tolMoles = 1e-18;
    propertiesY = {'T', 'p', 'h', 'e', 'g', 's', 'cp', 'cv', 'gamma_s', 'dVdT_p', 'dVdp_T', 'sound', 'W', 'uShock'};
    propertiesX = repmat({'equivalenceRatio'}, 1, length(propertiesY));
    
    % Define chemical system
    system = ChemicalSystem(database, listSpecies);
    
    % Initialize mixture
    mix = Mixture(system);
    
    % Define chemical state
    set(mix, {fuel}, 'fuel', 1);
    set(mix, {'N2', 'O2'}, 'oxidizer', [79, 21] / 21);
    
    % Define properties
    mixArray1 = setProperties(mix, 'temperature', 300, 'pressure', 1, 'equivalenceRatio',  value);
    
    % Initialize solver
    solver = DetonationSolver('problemType', 'DET', 'FLAG_RESULTS', false, 'tolMoles', tolMoles);
    
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