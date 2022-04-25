function problems_solved = run_validation_ROCKET_CEA_16
    % Run test validation_TP_CEA_1:
    % Contrasted with: NASA's Chemical Equilibrium with Applications software
    % Problem type: Equilibrium composition at exit of the rocket nozzle
    % Pressure chamber [bar] = 22;
    % Area ratio A_c/A_t = 2;
    % Area ratio A_e/A_t = 2.62;
    % Equivalence ratio [-] = 0.5:0.01:4
    % Initial mixture: RP-1 + LOX
    % List of species considered: ListSpecies('HC/O2/N2 PROPELLANTS')

    FLAG_SAVE = true;
    % Aratio = 1.1:0.1:2.62;
    Aratio = 2.62;
    for i = length(Aratio):-1:1
        % Inputs
        Fuel = 'RP_1';
        prefixDataName = 'RP1';
        filename = {strcat(prefixDataName, '_LOX_A16_ROCKET1.out'), strcat(prefixDataName, '_LOX_A16_ROCKET2.out')};
        LS =  'HC/O2/N2 PROPELLANTS';
        DisplaySpecies = {'CO2','CO','H2O','H2','O2','C2H2_acetylene',...
              'C2H4','C2H6','CH2CO_ketene','CH3','CH3CHO_ethanal','CH3OH',...
              'CH4','COOH','H','H2O2','HCHO_formaldehy','HCO','HCOOH','HO2',...
              'O','OH','Cbgrb'};
        tolN = 1e-18;
    
        % Load results CEA 
        results_CEA = data_CEA(filename, DisplaySpecies);
        % Combustion Toolbox
        results_CT = run_CT('ProblemType', 'ROCKET', 'TR', 298.15, 'pR', 22,...
                            'Species', LS,'S_Fuel', Fuel,'S_Oxidizer', 'O2bLb',...
                            'S_Inert', [], 'EquivalenceRatio', 1:0.01:4, 'tolN', tolN,...
                            'FLAG_IAC', false, 'Aratio_c', 2, 'Aratio', Aratio);
        problems_solved = length(results_CT.PD.range);
    
        Z(:, i) = cell2vector(results_CT.PS.strP, 'I_sp');
        OF(:, i) = cell2vector(results_CT.PS.strR, 'OF');
        
        % Display validation (plot)
        results_CT.C.mintol_display = 1e-10;
        results_CT.Misc.config.labelx = 'Mixture ratio $O/F$';
        results_CT.Misc.config.labely = 'Molar fraction $X_i$'; 
        results_CT.Misc.config.title = strcat('Area ratio $A_{\rm exit}/A_{\rm throat} = ', sprintf('%.2f', Aratio(i)), '$');
        % * Molar fractions
        fig1 = plot_molar_fractions_validation(results_CT, results_CEA, 'phi', 'Xi', DisplaySpecies, 'config', results_CT.Misc.config);
        set_title(gca, results_CT.Misc.config);
        % Save plots
        if FLAG_SAVE
            folderpath = strcat(pwd,'\Validations\Figures\');
            stack_trace = dbstack;
            filename = strcat('RP1_LOX_molar_', sprintf('%d', i));
            saveas(fig1, strcat(folderpath, filename), 'svg');
        end
    end
    
    % % Contour plot
    % [X, Y] = meshgrid(OF, Aratio);
    % ax = set_figure;
    % contourf(ax, X, Y, Z, 10);
    % c = colorbar;
    % xlabel(ax, 'Mixture ratio $O/F$');
    % ylabel(ax, 'Area ratio $A_{\rm exit}/A_{\rm throat}$');
    % c.Label.String = 'Specific impulse at sea level $I_{sp}$';
    % c.Label.Interpreter = 'latex';
    % set(ax, 'Layer','top')
end