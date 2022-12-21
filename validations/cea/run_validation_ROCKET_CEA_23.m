function problems_solved = run_validation_ROCKET_CEA_23
    % Run test validation_ROCKET_CEA_23:
    % Contrasted with: NASA's Chemical Equilibrium with Applications software
    % Problem type: Equilibrium composition at exit of the rocket nozzle
    % Pressure chamber [bar] = 101.325;
    % Model: Finite Area Chamber (FAC)
    % Area ratio A_c/A_t = 2;
    % Area ratio A_e/A_t = 8;
    % Equivalence ratio [-] = 0.5:0.01:4
    % Initial mixture: LH2 + LOX === H2bLb + O2bLb
    % List of species considered: list_species('HYDROGEN_L')

    % Inputs
    Fuel = 'CH6N2bLb';
    Aratio_c = 2;
    Aratio = 3;
    prefixDataName = 'MMH';
    filename = {[prefixDataName, '_N2O4_ROCKET_FAC1.out'], [prefixDataName, '_N2O4_ROCKET_FAC2.out']};
    LS =  'Soot formation extended';
    display_species = {'N2', 'H2O', 'O2', 'CO2', 'NO', 'OH', 'O', 'CO',...
        'H2', 'NO2', 'HO2', 'H', 'H2O2', 'N2O', 'HNO'};
    tolN = 1e-18;

    % Load results CEA 
    results_CEA = data_CEA(filename, display_species);
    % Combustion Toolbox
    results_CT = run_CT('ProblemType', 'ROCKET',...
                        'TR', 300,...
                        'pR', 101.325,...
                        'Species', LS,...
                        'S_Fuel', Fuel,...
                        'S_Oxidizer', 'N2O4',...
                        'ratio_oxidizers_O2', 1,...
                        'EquivalenceRatio', 0.5:0.01:4,...
                        'tolN', tolN,...
                        'FLAG_IAC', false,...
                        'Aratio_c', Aratio_c,...
                        'Aratio', Aratio);
    problems_solved = length(results_CT.PD.range);
    
    % Display validation (plot)
    results_CT.Misc.config.axis = 'tight';
    % * Molar fractions
    [~, fig1] = plot_molar_fractions(results_CT, results_CT.PS.strP, 'phi', 'Xi', 'validation', results_CEA, 'display_species', display_species);
    % * Properties mixture Exit - 1
    fig2 = plot_properties_validation(results_CT, results_CEA, {'phi', 'phi', 'phi', 'phi', 'phi', 'phi', 'phi', 'phi', 'phi'}, {'T', 'p', 'h', 'cP', 'cV', 'gamma_s', 'u', 'I_sp', 'I_vac'}, 'mix2');
    % Save plots
    folderpath = strcat(pwd,'\Validations\Figures\');
    stack_trace = dbstack;
    filename = stack_trace.name;
    saveas(fig1, strcat(folderpath, filename, '_molar'), 'svg');
    saveas(fig2, strcat(folderpath, filename, '_properties_1'), 'svg');
end