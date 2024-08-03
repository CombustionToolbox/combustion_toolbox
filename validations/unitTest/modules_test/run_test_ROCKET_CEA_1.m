function [max_rel_error_moles, max_rel_error_prop] = run_test_ROCKET_CEA_1(value, database)
    % Run test run_test_ROCKET_CEA_1:
    % Contrasted with: NASA's Chemical Equilibrium with Applications software
    % Problem type: Equilibrium composition at exit of the rocket nozzle
    % Pressure chamber [bar] = 101.325;
    % Model: Finite Area Chamber (FAC)
    % Area ratio A_c/A_t = 2;
    % Area ratio A_e/A_t = 3;
    % Equivalence ratio [-] = 0.5:0.01:4
    % Initial mixture: LCH4 + LOX === CH4bLb + O2bLb
    % List of species considered: list_species('Soot Formation Extended')
    
    % Import packages
    import combustiontoolbox.databases.NasaDatabase
    import combustiontoolbox.core.*
    import combustiontoolbox.rocket.*
    import combustiontoolbox.utils.display.*

    % Definitions
    fuel = 'CH4bLb';
    areaRatioChamber = 2;
    areaRatio = 3;
    prefixDataName = 'LCH4';
    filename = {strcat(prefixDataName, '_LOX_ROCKET_FAC1.out'), strcat(prefixDataName, '_LOX_ROCKET_FAC2.out')};
    listSpecies = 'Soot formation extended';
    displaySpecies = {'CO2','CO','H2O','H2','O2','C2H2_acetylene',...
          'C2H4','C2H6','CH2CO_ketene','CH3','CH3CHO_ethanal','CH3OH',...
          'CH4','COOH','H','H2O2','HCHO_formaldehy','HCO','HCOOH','HO2',...
          'O','OH','Cbgrb'};
    tolMoles = 1e-18;
    propertiesY = {'T', 'p', 'h', 'e', 'g', 's', 'cp', 'cv', 'gamma_s', 'dVdT_p', 'dVdp_T', 'sound', 'W', 'u', 'I_sp', 'I_vac'};
    propertiesX = repmat({'equivalenceRatio'}, 1, length(propertiesY));

    % Define chemical system
    system = ChemicalSystem(database, listSpecies);
    
    % Initialize mixture
    mix = Mixture(system);
    
    % Define chemical state
    set(mix, {fuel}, 'fuel', 1);
    set(mix, {'O2bLb'}, 'oxidizer', 1);
    
    % Define properties
    mixArray1 = setProperties(mix, 'temperature', 298.15, 'pressure', 101.325, 'equivalenceRatio', value, 'areaRatioChamber', areaRatioChamber, 'areaRatio', areaRatio);
    
    % Initialize solver
    solver = RocketSolver('problemType', 'ROCKET_FAC', 'tolMoles', tolMoles, 'FLAG_RESULTS', false);
    
    % Solve problem
    [~, ~, ~, ~, mixArray4] = solver.solveArray(mixArray1);

    % Load results CEA 
    resultsCEA = data_CEA(filename, displaySpecies);

    % Compute error: molar fractions
    max_rel_error_moles = compute_error_moles_CEA(mixArray4, resultsCEA, 'equivalenceRatio', value, 'Xi', displaySpecies, tolMoles);
    
    % Compute error: properties mixture 4
    max_rel_error_prop = compute_error_prop_CEA(mixArray4, resultsCEA, propertiesX, value, propertiesY, 'mix2');
end