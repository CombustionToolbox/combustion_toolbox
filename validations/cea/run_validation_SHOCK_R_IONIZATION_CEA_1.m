 function run_validation_SHOCK_R_IONIZATION_CEA_1
    % Run test validation_SHOCK_R_IONIZATION_CEA_1:
    % Contrasted with: NASA's Chemical Equilibrium with Applications software
    % Problem type: Planar reflected shock wave
    % Temperature [K]   = 300
    % Pressure    [bar] = 1
    % Incident velocity [m/s] = [~357, 8000]
    % Initial mixture: AIR (78.084% N2 + 20.9476% O2 + 0.9365% Ar + 0.0319% CO2)
    % List of species considered: list_species('Air_ions')
    
    % Import packages
    import combustiontoolbox.databases.NasaDatabase
    import combustiontoolbox.core.*
    import combustiontoolbox.shockdetonation.*
    import combustiontoolbox.utils.display.*

    % Benchmark?
    FLAG_BENCHMARK = false;

    % Definitions
    prefixDataName = 'airNASA_reflected_shocks_ionization';

    for i = 5:-1:1
        filename{i} = sprintf('%s%d.out', prefixDataName, i);
    end

    listSpecies = 'AIR_IONS';
    displaySpecies = {'eminus','Ar','Arplus','C','Cplus','Cminus','CN','CNplus','CNminus',...
                       'CNN','CO','COplus','CO2','CO2plus','C2','C2plus','C2minus',...
                       'CCN','CNC','OCCN','C2N2','C2O','C3','C3O2',...
                       'N','Nplus','Nminus','NCO','NO','NOplus','NO2','NO2minus',...
                       'NO3','NO3minus','N2','N2plus','N2minus','NCN','N2O','N2Oplus',...
                       'N2O3','N2O4','N2O5','N3','O','Oplus','Ominus','O2','O2plus',...
                       'O2minus','O3'};
    u1 = logspace(2, 5, 500); u1 = u1(u1<8000); u1 = u1(u1>=357);
    
    % Get Nasa database
    DB = NasaDatabase('FLAG_BENCHMARK', FLAG_BENCHMARK);
    
    % Define chemical system
    system = ChemicalSystem(DB, listSpecies);
    
    % Initialize mixture
    mix = Mixture(system);
    
    % Define chemical state
    set(mix, {'N2', 'O2', 'Ar', 'CO2'}, [78.084, 20.9476, 0.9365, 0.0319] / 20.9476);
    
    % Define properties
    mixArray1 = setProperties(mix, 'temperature', 300, 'pressure', 1, 'u1', u1);
    
    % Initialize solver
    solver = ShockSolver('problemType', 'SHOCK_R', 'FLAG_RESULTS', false);
    
    % Solve problem
    [mixArray1, mixArray2, mixArray3] = solver.solveArray(mixArray1);
    
    if FLAG_BENCHMARK
        return
    end
    
    % Prepare data
    for i = 1:length(mixArray3)
        mixArray3(i).W = mixArray3(i).W * 1e3; % [g/mol]
    end

    % Load results CEA 
    resultsCEA = data_CEA(filename, displaySpecies);
    
    % Plot molar fractions
    fig1 = plotComposition(mixArray3(1), mixArray1, 'u', 'Xi', 'mintol', 1e-5, 'y_var', mixArray3);

    % Properties mixture 1
    fig2 = plotProperties(repmat({'u_preshock'}, 1, 8), [mixArray1.u], {'T', 'p', 'rho', 'h', 'e', 'g', 's', 'gamma_s'}, mixArray1, 'basis', {[], [], [], 'mi', 'mi', 'mi', 'mi', []}, 'validation', resultsCEA.mix1);

    % Properties mixture 2 - 1
    fig3 = plotProperties(repmat({'u_preshock'}, 1, 8), [mixArray1.u], {'T', 'p', 'rho', 'h', 'e', 'g', 's', 'gamma_s'}, mixArray3, 'basis', {[], [], [], 'mi', 'mi', 'mi', 'mi', []}, 'validation', resultsCEA.mix2);

    % Properties mixture 2 - 2
    fig4 = plotProperties(repmat({'u_preshock'}, 1, 6), [mixArray1.u], {'cp', 'cv', 'dVdT_p', 'dVdp_T', 'sound', 'W'}, mixArray3, 'basis', {'mi', 'mi', [], [], [], []}, 'validation', resultsCEA.mix2);

    % Save plots
    folderpath = fullfile(pwd, 'validations', 'figures');
    stack_trace = dbstack;
    filename = stack_trace.name;
    saveas(fig1, fullfile(folderpath, strcat(filename, '_molar')), 'svg');
    saveas(fig2, fullfile(folderpath, strcat(filename, '_properties_1')), 'svg');
    saveas(fig3, fullfile(folderpath, strcat(filename, '_properties_2')), 'svg');
    saveas(fig4, fullfile(folderpath, strcat(filename, '_properties_3')), 'svg');
 end