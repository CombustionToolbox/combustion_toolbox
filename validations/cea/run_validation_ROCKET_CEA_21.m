function problems_solved = run_validation_ROCKET_CEA_21
    % Run test validation_ROCKET_CEA_21:
    % Contrasted with: NASA's Chemical Equilibrium with Applications software
    % Problem type: Equilibrium composition at exit of the rocket nozzle
    % Pressure chamber [bar] = 101.325;
    % Model: Finite Area Chamber (FAC)
    % Area ratio A_c/A_t = 2;
    % Area ratio A_e/A_t = 3;
    % Equivalence ratio [-] = 0.5:0.01:4
    % Initial mixture: LCH4 + LOX === CH4bLb + O2bLb
    % List of species considered: list_species('Soot Formation Extended')

    % Inputs
    Fuel = 'CH4bLb';
    Aratio_c = 2;
    Aratio = 3;
    prefixDataName = 'LCH4';
    filename = {[prefixDataName, '_LOX_ROCKET_FAC1.out'], [prefixDataName, '_LOX_ROCKET_FAC2.out']};
    LS =  'Soot Formation Extended';
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
                        'pR', 101.325,...
                        'Species', LS,...
                        'S_Fuel', Fuel,...
                        'S_Oxidizer', 'O2bLb',...
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