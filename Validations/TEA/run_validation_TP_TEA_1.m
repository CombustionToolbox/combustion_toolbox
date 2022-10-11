function problems_solved = run_validation_TP_TEA_1
    % VALIDATION: TP_TEA
    %
    % Compute equilibrium composition at defined temperature and pressure.
    % Reproduce the example case of TEA by Jasmina Blecic.
    % URL RESULTS TEA:
    % https://github.com/dzesmin/TEA/tree/master/doc/examples/quick_example/results 
    %   
    %
    % @author: Alberto Cuadra Lara
    %          PhD Candidate - Group Fluid Mechanics
    %          Universidad Carlos III de Madrid
    %                  
    % Last update April 15 2022
    
    % Inputs
    load Validation_TP_TEA_1 Pressure Temp results_TEA

    metallicity = 1;
    Fuel = {'H', 'He', 'C', 'N', 'O'};
    N_Fuel = abundances2moles(Fuel, 'abundances.txt', metallicity);

    LS = {'C', 'CH4', 'CO2', 'CO', 'H2', 'H', 'H2O', 'He', 'N2', 'N', 'NH3', 'O'};
    T = linspace(Temp(1), Temp(end), 300);
    p = logspace(-5, 2, 300);
    Oxidizer = {};
    % Tunning paramenters
    tolN = 1e-32;
    % Custom Plots 
    DisplaySpecies = LS;
    % Combustion Toolbox
    results_CT = run_CT('ProblemType', 'TP',...
                        'Temperature', T,...
                        'Pressure', p, ...
                        'Species', LS,...
                        'S_Fuel', Fuel,...
                        'N_Fuel', N_Fuel,...
                        'S_Oxidizer', Oxidizer,...
                        'tolN', tolN);
    problems_solved = length(results_CT.PD.range);
    % Display validation (plot)
    % * Molar fractions
    fig1 = plot_molar_fractions_validation(results_CT, results_TEA, 'T', 'Xi', DisplaySpecies, 'tolN', tolN, 'nfrec', 3);
    % Save plots
    folderpath = strcat(pwd,'\Validations\Figures\');
    stack_trace = dbstack;
    filename = stack_trace.name;
    saveas(fig1, strcat(folderpath, filename, '_molar'), 'svg');
end
