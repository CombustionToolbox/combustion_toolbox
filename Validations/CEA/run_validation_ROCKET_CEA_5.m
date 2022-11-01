function problems_solved = run_validation_ROCKET_CEA_5
    % Run test validation_ROCKET_CEA_5:
    % Contrasted with: NASA's Chemical Equilibrium with Applications software
    % Problem type: Equilibrium composition at exit of the rocket nozzle
    % Pressure chamber [bar] = 22;
    % Area ratio A_c/A_t = 2;
    % Area ratio A_e/A_t = 1.5;
    % Equivalence ratio [-] = 0.5:0.01:4
    % Initial mixture: RP-1 + LOX
    % List of species considered: list_species('HC/O2/N2 PROPELLANTS')

    % Inputs
    Fuel = 'RP_1';
    Aratio = 1.5;
    prefixDataName = 'RP1';
    filename = {strcat(prefixDataName, '_LOX_A5_ROCKET1.out'), strcat(prefixDataName, '_LOX_A5_ROCKET2.out')};
    LS =  'HC/O2/N2 PROPELLANTS';
    display_species = {'CO2','CO','H2O','H2','O2','C2H2_acetylene',...
          'C2H4','C2H6','CH2CO_ketene','CH3','CH3CHO_ethanal','CH3OH',...
          'CH4','COOH','H','H2O2','HCHO_formaldehy','HCO','HCOOH','HO2',...
          'O','OH','Cbgrb'};
    tolN = 1e-18;

    % Load results CEA 
    results_CEA = data_CEA(filename, display_species);
    % Combustion Toolbox
    results_CT = run_CT('ProblemType', 'ROCKET',...
                        'TR', 298.15,...
                        'pR', 22,...
                        'Species', LS,...
                        'S_Fuel', Fuel,...
                        'S_Oxidizer', 'O2bLb',...
                        'ratio_oxidizers_O2', 1,...
                        'EquivalenceRatio', 1:0.01:4,...
                        'tolN', tolN,...
                        'FLAG_IAC', false,...
                        'Aratio_c', 2,...
                        'Aratio', Aratio);
    problems_solved = length(results_CT.PD.range);
    
    % Display validation (plot)
    results_CT.Misc.config.labelx = 'Mixture ratio $O/F$';
    results_CT.Misc.config.labely = 'Molar fraction $X_i$'; 
    results_CT.Misc.config.title = strcat('Area ratio $A_{\rm exit}/A_{\rm throat} = ', sprintf('%.2f', Aratio), '$');
    % * Molar fractions
%         results_CEA.OF = OF(:, i);
    fig1 = plot_molar_fractions_validation(results_CT, results_CEA, 'phi', 'Xi', display_species, 'config', results_CT.Misc.config);
    set_title(gca, results_CT.Misc.config);
    % Save plots
    folderpath = strcat(pwd,'\Validations\Figures\');
    stack_trace = dbstack;
    filename = stack_trace.name;
    saveas(fig1, strcat(folderpath, filename, '_molar'), 'svg');
end