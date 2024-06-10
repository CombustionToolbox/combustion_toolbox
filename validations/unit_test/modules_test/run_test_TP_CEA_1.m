function [max_rel_error_moles, max_rel_error_prop] = run_test_TP_CEA_1(value, database)
    % Run test_TP_CEA_1:
    % Contrasted with: NASA's Chemical Equilibrium with Applications software
    % Problem type: Equilibrium composition at defined T and p
    % Temperature [K]   = value;
    % Pressure    [bar] = 1.01325;
    % Initial mixture: Si + 9 C6H5OH_phenol
    % List of species considered: All (see routine find_products)
    
    % Import packages
    import combustiontoolbox.databases.NasaDatabase
    import combustiontoolbox.core.*
    import combustiontoolbox.equilibrium.*
    import combustiontoolbox.utils.display.*

    % Definitions
    species = {'Si', 'C6H5OH_phenol'};
    moles = [1, 9];
    prefixDataName = 'C6H5OH_phenol_and_Si';
    filename = {strcat(prefixDataName, '_TP1.out')};
    displaySpecies = {'H2', 'H', 'CH4', 'C2H2_acetylene', 'SiC2',...
        'Si', 'C', 'C2H', 'C3', 'CO', 'CO2', 'Cbgrb', 'SiCbbb', 'H2O',...
        'SiO2ba_qzb','SiO2bb_qzb','SiO2bb_crtb','H2ObLb','H2Obcrb'};
    tolMoles = 1e-18;
    propertiesY = {'p', 'h', 'e', 'g', 's', 'cp', 'cv', 'gamma_s', 'dVdT_p', 'dVdp_T', 'sound', 'W'};
    propertiesX = repmat({'T'}, 1, length(propertiesY));

    % Define chemical system
    system = ChemicalSystem(database);
    
    % Initialize mixture
    mix = Mixture(system);
    
    % Define chemical state
    set(mix, species, moles);
    
    % Define properties
    mixArray = setProperties(mix, 'temperature', value, 'pressure', 1 * 1.01325);
    
    % Initialize solver
    solver = EquilibriumSolver('problemType', 'TP', 'FLAG_RESULTS', false, 'tolMoles', tolMoles);
    
    % Solve problem
    solver.solveArray(mixArray);
    
    % Load results CEA 
    results_CEA = data_CEA(filename, displaySpecies);

    % Compute error: molar fractions
    max_rel_error_moles = compute_error_moles_CEA(mixArray, results_CEA, 'T', value, 'Xi', displaySpecies, tolMoles);
    
    % Compute error: properties mixture 2
    max_rel_error_prop = compute_error_prop_CEA(mixArray, results_CEA, propertiesX, value, propertiesY, 'mix2');
end