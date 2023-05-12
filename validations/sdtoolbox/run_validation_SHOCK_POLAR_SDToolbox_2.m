function problems_solved = run_validation_SHOCK_POLAR_SDToolbox_2
    % Run test validation_SHOCK_POLAR_SDToolbox_2
    % Contrasted with: Caltech's SD Toolbox and CANTERA
    % Problem type: Shock polar
    % Temperature [K]   = 300
    % Pressure    [bar] = 1.01325
    % Initial mixture: AIR (79% N2 + 21% O2)
    % List of species considered: list_species('AIR_IONS')
    
    % Inputs

%     Oxidizer = {'N2', 'O2', 'Ar', 'CO2'};
%     moles = [3.7276, 1.0000, 0.0447, 0.0015];

    Oxidizer = {'N2', 'O2', 'Ar'};
    moles = [78, 21, 1]/21;
%     LS = 'AIR_IONS';

    LS = {'O2', 'O', 'N2', 'N', 'NO', 'eminus', 'Nplus', 'NOplus',...
          'N2plus', 'Oplus', 'O2plus', 'Ar'};
    
    u1 = 3.472107491008314e+02 * [2, 3, 5, 14];
    filename = 'shock_polar_equilwithions_ar_SDToolbox';
    % Tunning parameters
    tolN = 1e-14;
    N_points_polar = 300;
    % Combustion Toolbox
    results_CT = run_CT('ProblemType', 'SHOCK_POLAR',...
                        'Species', LS,...
                        'S_Oxidizer', Oxidizer,...
                        'N_Oxidizer', moles,...
                        'u1', u1,...
                        'tolN', tolN,...
                        'N_points_polar', N_points_polar);
    problems_solved = get_problems_solved(results_CT.PS.strP, 'polar', 'theta');
    % Load results SDToolbox 
    results_SDToolbox = load_struct(filename, 'data');
    % Display validation (plot)
    [fig1, fig2] = plot_validation_shock_polar_SDToolbox(results_CT, results_SDToolbox, results_CT.Misc.config);
    % Save plots
    folderpath = strcat(pwd,'\Validations\Figures\');
    stack_trace = dbstack;
    filename = stack_trace.name;
    saveas(fig1, strcat(folderpath, filename, '_properties_1'), 'svg');
    saveas(fig2, strcat(folderpath, filename, '_properties_2'), 'svg');
end