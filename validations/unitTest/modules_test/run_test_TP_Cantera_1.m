function [max_rel_error_moles, max_rel_error_prop] = run_test_TP_Cantera_1(value, database)
    % Run test_TP_CEA_1:
    % Contrasted with: Cantera
    % Problem type: Equilibrium composition at defined T and p
    % Temperature [K]   = 200;
    % Pressure    [bar] = value;
    % Initial mixture: AIR_IDEAL (79% N2 + 21% O2)
    % Equation of state: Peng-Robinson
    
    % Import packages
    import combustiontoolbox.databases.NasaDatabase
    import combustiontoolbox.core.*
    import combustiontoolbox.equilibrium.*
    import combustiontoolbox.utils.display.*

    % Definitions
    species = {'N2', 'O2'};
    moles = [79/21, 1];
    prefixDataName = 'air';
    filename = strcat('Cantera_', prefixDataName, '_TP1.txt');
    propertiesY = {'rho', 'h', 'e', 'g', 's', 'cp', 'cv', 'gamma_s', 'dVdp_T', 'dVdT_p', 'dPdV_T', 'dPdT_V'};
    propertiesX = repmat({'p'}, 1, length(propertiesY));

    % Add Peng-Robinson parameters
    database.species = addPengRobinsonProperties(database.species);

    % Define chemical system
    system = ChemicalSystem(database);
    
    % Define equation of state
    eos = EquationStatePengRobinson();

    % Initialize mixture
    mix = Mixture(system, 'eos', eos);
    
    % Define chemical state
    set(mix, species, moles);
    
    % Define properties
    mixArray = setProperties(mix, 'temperature', 200, 'pressure', value);
    
    % Load results Cantera
    results_Cantera = readtable(filename, 'PreserveVariableNames', true);

    % Compute error: molar fractions
    max_rel_error_moles = 0;
    
    % Compute error: properties mixture
    max_rel_error_prop = compute_error_prop_Cantera(mixArray, results_Cantera, propertiesX, value, propertiesY);
end