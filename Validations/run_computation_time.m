function [computational_time, problems_solved] = run_computation_time(ProblemType, S_Fuel, LS_0, N_LS_min, N_LS_frec, N_repetitions)
    % Run test validation_DET_CEA_1:
    % Contrasted with: NASA's Chemical Equilibrium with Applications software
    % Problem type: Chapman-Jouguet Detonation
    % Temperature [K]   = 300;
    % Pressure    [bar] = 1;
    % Equivalence ratio [-] = 0.5:0.01:4
    % Initial mixture: C2H2_acetylene + AIR_IDEAL (79% N2 + 21% O2)
    % List of species considered: ListSpecies('Soot Formation Extended')
    % 
    % Args:
    %
    %
    % Returns:



    % Default values
    LS_0 
    N_LS_min = 9; % Minimum number of species to consider
    N_LS_frec = 1; % Step for chemical species
    N_repetitions = 3; % Repetitions of the routine

    




    % Inputs 
%     LS =  'Soot Formation Extended';
%     DisplaySpecies = {'CO2', 'CO', 'H2O', 'H2', 'O2', 'N2',...
%                       'HCN','H','OH','O','CN','NH3','CH4','C2H4','CH3',...
%                       'NO','HCO','NH2','NH','N','CH','Cbgrb'};
    tolN = 1e-18;
    
    N_LS_max = length(LS_0);
    LS_min = LS_0(1:N_LS_min);
    Misc = Miscellaneous();
    config = Misc.config;
    config.labelx = 'Number of species';
    config.labely = 'Computational time [s]';
    ax = set_figure(config);
    for j = N_repetitions:-1:1
        for i = N_LS_max:-N_LS_frec:N_LS_min + 1
            LS = [LS_min, LS_0(N_LS_min+1:N_LS_frec:end)];
        
            % Computations
            results_CT = run_CT('ProblemType', ProblemType,...
                                'Species', LS,...
                                'S_Fuel', S_Fuel,...
                                'EquivalenceRatio', 0.5:0.5:5,...
                                'tolN', tolN);
            computational_time(i) = results_CT.Misc.timer_loop;
        end
        problems_solved = length(N_LS_max:-N_LS_frec:N_LS_min + 1) * length(results_CT.PD.range);
        
        % Set dataset (number of species, CPU time in seconds)
        x = N_LS_max:-N_LS_frec:N_LS_min + 1;
        y(j, :) = computational_time(N_LS_min+1:end);
        % Plot scatter data and linear regression
        scatter(ax, x, y(j, :), N_LS_max:-N_LS_frec:N_LS_min + 1,'MarkerEdgeColor',[0 .5 .5],...
                  'MarkerFaceColor',[0 .7 .7],...
                  'LineWidth',1.5)
    end
    y_mean = mean(y);
    % Obtain polynomial regression of order 1
    y_poly = polynomial_regression(x, y_mean, 1);
    % Plot polynomial regression
    plot(ax, x, y_poly, 'LineWidth', results_CT.Misc.config.linewidth);
end