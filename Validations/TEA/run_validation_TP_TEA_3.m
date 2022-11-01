function problems_solved = run_validation_TP_TEA_3
    % VALIDATION: TP_TEA_3
    %
    % Compute equilibrium composition at defined temperature and pressure.
    % Reproduce Fig.6 (b) "Thermochemical equilibrium vertical distributions
    % with a metallicity 10 of WASP-43b assuming the T-P profile in Fig 4"
    % [TEA by Jasmina Blecic].
    %
    % URL RESULTS TEA:
    % https://github.com/dzesmin/RRC-BlecicEtal-2015a-ApJS-TEA/tree/master/Fig6/WASP43b-10xsolar
    %   
    %
    % @author: Alberto Cuadra Lara
    %          PhD Candidate - Group Fluid Mechanics
    %          Universidad Carlos III de Madrid
    %                  
    % Last update Oct 12 2022
    
    % Inputs
    load Validation_TP_TEA_3 Pressure Temp results_TEA

    metallicity = 10;
    Fuel = {'H', 'He', 'C', 'N', 'O', 'S'};
    N_Fuel = abundances2moles(Fuel, 'abundances.txt', metallicity);
    Oxidizer = {};
    LS = {'C2H2_acetylene', 'C2H4', 'C', 'CH4', 'CO2', 'CO', 'H2',...
          'H2O', 'H2S', 'H', 'HCN', 'He', 'HS_M', 'N2', 'N', 'NH3',...
          'O', 'S'};
    
    T = Temp;
    p = Pressure;
    % Tunning paramenters
    tolN = 1e-32;
    % Custom Plots 
    display_species = LS;
    mintol = 1e-21;
    % Combustion Toolbox
    results_CT = run_CT('ProblemType', 'TP',...
                        'Temperature', T,...
                        'Pressure', p,...
                        'Species', LS,...
                        'S_Fuel', Fuel,...
                        'N_Fuel', N_Fuel,...
                        'S_Oxidizer', Oxidizer,...
                        'tolN', tolN);
    problems_solved = length(results_CT.PD.range);
    % Display validation (plot)
    % * Molar fractions
    [~, fig1] = plot_molar_fractions(results_CT, results_CT.PS.strP, ...
        'Xi', 'p', 'validation', results_TEA, 'nfrec', 3,...
        'ydir', 'reverse', 'xscale', 'log', 'mintol', mintol,...
        'display_species', display_species);
    % Save plots
    folderpath = strcat(pwd,'\Validations\Figures\');
    stack_trace = dbstack;
    filename = stack_trace.name;
    saveas(fig1, strcat(folderpath, filename, '_molar'), 'svg');
end
