function problems_solved = run_validation_TP_TEA_4
    % VALIDATION: TP_TEA_3
    %
    % Compute equilibrium composition at defined temperature and pressure.
    % Reproduce Fig.6 (b) "Thermochemical equilibrium vertical distributions
    % with a metallicity 50 of WASP-43b assuming the T-P profile in Fig 4"
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
    % Last update April 15 2022
    
    % Inputs
    load Validation_TP_TEA_4 Pressure Temp results_TEA
    Fuel = {'H', 'He', 'C', 'N', 'O', 'S'};
    Xi_abundances = abundances2moles(Fuel, 'abundances_WASP43b_50xsolar.txt')';
    N_Fuel = Xi_abundances;
    Oxidizer = {}; Inert = {};
    LS = {'C2H2_acetylene', 'C2H4', 'C', 'CH4', 'CO2', 'CO', 'H2', 'H2O', 'H2S', 'H', 'HCN', 'He', 'HS_M', 'N2', 'N', 'NH3', 'O', 'S'};
    
    T = Temp;
    p = Pressure;
    % Tunning paramenters
    tolN = 1e-30;
    % Custom Plots 
    DisplaySpecies = LS;
    mintol = 1e-21;
    % Combustion Toolbox
    results_CT = run_CT('ProblemType', 'TP', 'Temperature', T, 'Pressure', p, ...
        'Species', LS, 'S_Fuel', Fuel, 'N_Fuel', N_Fuel, 'S_Oxidizer', Oxidizer,...
        'S_Inert', Inert, 'tolN', tolN);
    problems_solved = length(results_CT.PD.range);
    % Display validation (plot)
    % * Molar fractions
    fig1 = plot_molar_fractions_validation(results_CT, results_TEA, 'Xi', 'p', DisplaySpecies, 'mintol', mintol, 'nfrec', 3,...
        'ydir', 'reverse', 'xscale', 'log');
    % Save plots
    folderpath = strcat(pwd,'\Validations\Figures\');
    stack_trace = dbstack;
    filename = stack_trace.name;
    saveas(fig1, strcat(folderpath, filename, '_molar'), 'svg');

%     results_TEA.Xi = [C2H2_g, C2H4_g, C_g, CH4_g, CO2_g, CO_g, H2_ref, H2O_g, H2S_g, H_g, HCN_g, He_ref, HS_g, N2_ref, N_g, NH3_g, O_g, S_g]';
%     results_TEA.T = Temp;
%     results_TEA.p = Pressure;
end
