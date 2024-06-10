 function [max_rel_error_prop_mix1, max_rel_error_prop_mix3] = run_test_SHOCK_R_IONIZATION_CEA_1(value, database)
    % Run test_SHOCK_IONIZATION_CEA_1:
    % Contrasted with: NASA's Chemical Equilibrium with Applications software
    % Problem type: Planar incident shock wave
    % Temperature [K]   = 300
    % Pressure    [bar] = 1
    % Incident velocity [m/s] = value
    % Initial mixture: AIR (78.084% N2 + 20.9476% O2 + 0.9365% Ar + 0.0319% CO2)
    % List of species considered: All (see method findProducts from ChemicalSystem class)
    
    % Import packages
    import combustiontoolbox.databases.NasaDatabase
    import combustiontoolbox.core.*
    import combustiontoolbox.shockdetonation.*
    import combustiontoolbox.utils.display.*

    % Definitions
    prefixDataName = 'airNASA_reflected_shocks_ionization';

    for i = 5:-1:1
        filename{i} = sprintf('%s%d.out', prefixDataName, i);
    end

    displaySpecies = {'eminus','Ar','Arplus','C','Cplus','Cminus','CN','CNplus','CNminus',...
                      'CNN','CO','COplus','CO2','CO2plus',...
                      'N','Nplus','Nminus','NCO','NO','NOplus','NO2','NO2minus',...
                      'N2','N2plus','N2minus','NCN','N2O','N2Oplus','N3',...
                      'O','Oplus','Ominus','O2','O2plus','O2minus','O3'};
    
    % Define chemical system
    system = ChemicalSystem(database, [], 'FLAG_ION', true);
    
    % Initialize mixture
    mix = Mixture(system);
    
    % Define chemical state
    set(mix, {'N2', 'O2', 'Ar', 'CO2'}, [78.084, 20.9476, 0.9365, 0.0319] / 20.9476);
    
    % Define properties
    mixArray1 = setProperties(mix, 'temperature', 300, 'pressure', 1, 'u1', value);
    
    % Initialize solver
    solver = ShockSolver('problemType', 'SHOCK_R', 'FLAG_RESULTS', false);
    
    % Solve problem
    [mixArray1, mixArray2, mixArray3] = solver.solveArray(mixArray1);

    % Load results CEA 
    resultsCEA = data_CEA(filename, displaySpecies);

    % Compute error: molar fractions
    % max_rel_error_moles = compute_error_moles_CEA(results_CT, results_CEA, 'u', value, 'Xi', DisplaySpecies);
    
    % Compute error: Properties mixture 1
    propertiesY = {'T', 'p', 'h', 'e', 'g', 's', 'cp', 'cv', 'gamma_s', 'sound', 'W'};
    propertiesX = repmat({'u'}, 1, length(propertiesY));
    max_rel_error_prop_mix1 = compute_error_prop_CEA(mixArray1, resultsCEA, propertiesX, value, propertiesY, 'mix1');
    
    % Compute error: Properties mixture 3
    propertiesY = {'T', 'p', 'h', 'e', 'g', 's', 'cp', 'cv', 'gamma_s', 'dVdT_p', 'dVdp_T', 'sound', 'W'};
    propertiesX = repmat({'u'}, 1, length(propertiesY));
    max_rel_error_prop_mix3 = compute_error_prop_CEA(mixArray3, resultsCEA, propertiesX, value, propertiesY, 'mix2');
end