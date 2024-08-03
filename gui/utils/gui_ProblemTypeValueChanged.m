function gui_ProblemTypeValueChanged(app)
    % Clear GUI results tab (except UITree) and update GUI items for the problem selected

    % 1. Clear results tab (except UITree)
    gui_clear_results(app);

    % 2. Update GUI items depending of the problem selected    
    switch app.ProblemType.Value
        case 'TP' % * TP: Equilibrium composition at defined T and p
            % Visible flags
            app.FLAG_IAC.Visible = 'off';
            app.FLAG_IAC.Value = true;
            gui_FLAG_IACValueChanged(app);
            % Update input items
            app.PP1.Visible = 'on'; app.PP2.Visible = 'on';
            app.text_P1.Visible = 'on';
            app.text_RP2.Text = 'Pressure [bar]';
            % Update Additional constraints panel
            app.AdditionalconstraintsPanel.Visible = 'off';
            app.AdditionalconstraintsPanel.Title = 'Additional constraints';
            % Set default input values
            app.PR1.Value = '300'; 
            app.PR2.Value = '1';
            app.PP1.Value = '2500';
            app.PP2.Value = app.PR2.Value;
            % Set invisible shock/detonation items
            gui_visible_shocks(app, false);
            % Set invisible oblique shock/detontions items
            gui_visible_oblique(app, false);
            % Set invisible rocket items
            gui_visible_rocket(app, false);
        case 'HP' % * HP: Adiabatic T and composition at constant p
            % Visible flags
            app.FLAG_IAC.Visible = 'off';
            app.FLAG_IAC.Value = true;
            gui_FLAG_IACValueChanged(app);
            % Update input items
            app.PP1.Visible = 'off'; app.PP2.Visible = 'on';
            app.PR3.Visible = 'off'; app.PR4.Visible = 'off';
            app.PP3.Visible = 'off'; app.PP4.Visible = 'off';
            app.PR5.Visible = 'off'; app.PP5.Visible = 'off';
            app.text_P1.Visible = 'on';
            app.text_RP2.Text = 'Pressure [bar]';
            % Visible flags
            app.FLAG_IAC.Visible = 'off';
            % Update Additional constraints panel
            app.AdditionalconstraintsPanel.Visible = 'on';
            app.AdditionalconstraintsPanel.Title = 'Additional constraints';
            app.text_RP3.Text = 'Constant Enthalpy: hP = hR';
            app.text_RP.Visible = 'on'; app.text_RP.Text = 'Products';
            app.text_R2.Visible = 'off'; app.text_P2.Visible = 'off';
            app.text_R2.Text = 'Reactants'; app.text_P2.Text = 'Products';
            app.text_RP4.Visible = 'off';
            app.text_RP5.Visible = 'off';
            % Set default input values
            app.PR1.Value = '300';
            app.PR2.Value = '1';
            app.PP2.Value = app.PR2.Value;
            % Set invisible shock/detonation items
            gui_visible_shocks(app, false);
            % Set invisible oblique shock/detontions items
            gui_visible_oblique(app, false);
            % Set invisible rocket items
            gui_visible_rocket(app, false);
        case 'SP' % * SP: Isentropic (i.e., adiabatic) compression/expansion to a specified p
            % Visible flags
            app.FLAG_IAC.Visible = 'off';
            app.FLAG_IAC.Value = true;
            gui_FLAG_IACValueChanged(app);
            % Update input items
            app.PP1.Visible = 'off'; app.PP2.Visible = 'on';
            app.PR3.Visible = 'off'; app.PR4.Visible = 'off';
            app.PP3.Visible = 'off'; app.PP4.Visible = 'off';
            app.PR5.Visible = 'off'; app.PP5.Visible = 'off';
            app.text_P1.Visible = 'on';
            app.text_RP2.Text = 'Pressure [bar]';
            % Visible flags
            app.FLAG_IAC.Visible = 'off';
            % Update Additional constraints panel
            app.AdditionalconstraintsPanel.Visible = 'on';
            app.AdditionalconstraintsPanel.Title = 'Additional constraints';
            app.text_RP.Visible = 'on'; app.text_RP.Text = 'Products';
            app.text_R2.Visible = 'off'; app.text_P2.Visible = 'off';
            app.text_R2.Text = 'Reactants'; app.text_P2.Text = 'Products';
            app.text_RP3.Text = 'Constant Entropy: SP = SR';
            app.text_RP4.Visible = 'off';
            app.text_RP5.Visible = 'off';
            % Set default input values
            app.PR1.Value = '300';
            app.PR2.Value = '1';
            app.PP2.Value = '100';
            % Set invisible shock/detonation items
            gui_visible_shocks(app, false);
            % Set invisible oblique shock/detontions items
            gui_visible_oblique(app, false);
            % Set invisible rocket items
            gui_visible_rocket(app, false);
        case 'TV' % * TV: Equilibrium composition at defined T and constant v
            % Visible flags
            app.FLAG_IAC.Visible = 'off';
            app.FLAG_IAC.Value = true;
            gui_FLAG_IACValueChanged(app);
            % Update input items
            app.PP1.Visible = 'on'; app.PP2.Visible = 'on'; 
            app.PR3.Visible = 'off'; app.PR4.Visible = 'off';
            app.PP3.Visible = 'off'; app.PP4.Visible = 'off';
            app.PR5.Visible = 'off'; app.PP5.Visible = 'off';
            app.text_P1.Visible = 'on';
            app.text_RP2.Text = 'Specific volume [m3/kg]';
            % Visible flags
            app.FLAG_IAC.Visible = 'off';
            % Update Additional constraints panel
            app.AdditionalconstraintsPanel.Visible = 'off';
            app.AdditionalconstraintsPanel.Title = 'Additional constraints';
            % Set default input values
            app.PR1.Value = '300';
            app.PR2.Value = '1';
            app.PP1.Value = '2500';
            app.PP2.Value = app.PR2.Value;
            % Set invisible shock/detonation items
            gui_visible_shocks(app, false);
            % Set invisible oblique shock/detontions items
            gui_visible_oblique(app, false);
            % Set invisible rocket items
            gui_visible_rocket(app, false);
        case 'EV' % * EV: Equilibrium composition at Adiabatic T and constant v
            % Visible flags
            app.FLAG_IAC.Visible = 'off';
            app.FLAG_IAC.Value = true;
            gui_FLAG_IACValueChanged(app);
            % Update input items
            app.PP1.Visible = 'off'; app.PP2.Visible = 'off'; 
            app.PR3.Visible = 'off'; app.PR4.Visible = 'off';
            app.PP3.Visible = 'off'; app.PP4.Visible = 'off';
            app.PR5.Visible = 'off'; app.PP5.Visible = 'off'; 
            app.text_P1.Visible = 'on';
            app.text_RP2.Text = 'Specific volume [m3/kg]';
            % Visible flags
            app.FLAG_IAC.Visible = 'off';
            % Update Additional constraints panel
            app.AdditionalconstraintsPanel.Visible = 'on';
            app.AdditionalconstraintsPanel.Title = 'Additional constraints';
            app.text_RP.Visible = 'on'; app.text_RP.Text = 'Products';
            app.text_R2.Visible = 'off'; app.text_P2.Visible = 'off';
            app.text_R2.Text = 'Reactants'; app.text_P2.Text = 'Products';
            app.text_RP3.Text = 'Constant Internal energy: eP = eR';
            app.text_RP4.Visible = 'on'; app.text_RP4.Text = 'Constant Volume: vP = vR';
            app.text_RP5.Visible = 'off';
            % Set default input values
            app.PR1.Value = '300';
            app.PR2.Value = '1';
            app.PP2.Value = app.PR2.Value;
            % Set invisible shock/detonation items
            gui_visible_shocks(app, false);
            % Set invisible oblique shock/detontions items
            gui_visible_oblique(app, false);
            % Set invisible rocket items
            gui_visible_rocket(app, false);
        case 'SV' % * SV: Isentropic (i.e., fast adiabatic) compression/expansion to a specified v
            % Visible flags
            app.FLAG_IAC.Visible = 'off';
            app.FLAG_IAC.Value = true;
            gui_FLAG_IACValueChanged(app);
            % Update input items
            app.PP1.Visible = 'off'; app.PP2.Visible = 'off'; 
            app.PR3.Visible = 'off'; app.PR4.Visible = 'off';
            app.PP3.Visible = 'off'; app.PP4.Visible = 'on';
            app.PR5.Visible = 'off'; app.PP5.Visible = 'off';
            app.text_P1.Visible = 'on';
            app.text_RP2.Text = 'Specific volume [m3/kg]';
            % Visible flags
            app.FLAG_IAC.Visible = 'off';
            % Update Additional constraints panel
            app.AdditionalconstraintsPanel.Visible = 'on';
            app.AdditionalconstraintsPanel.Title = 'Additional constraints';
            app.text_RP.Visible = 'on'; app.text_RP.Text = 'Products';
            app.text_R2.Visible = 'off'; app.text_P2.Visible = 'off';
            app.text_R2.Text = 'Reactants'; app.text_P2.Text = 'Products';
            app.text_RP3.Text = 'Constant Entropy: SP = SR';
            app.text_RP4.Visible = 'on'; app.text_RP4.Text = 'Volume Products/Reactants';
            app.text_RP5.Visible = 'off';
            % Set default input values
            app.PR1.Value = '300';
            app.PR2.Value = '1';
            app.PP2.Value = app.PR2.Value;
            app.PP4.Value = '2';
            % Set invisible shock/detonation items
            gui_visible_shocks(app, false);
            % Set invisible oblique shock/detontions items
            gui_visible_oblique(app, false);
            % Set invisible rocket items
            gui_visible_rocket(app, false);
        case 'SHOCK_I' % * SHOCK_I: CALCULATE PLANAR INCIDENT SHOCK WAVE
            % Visible flags
            app.FLAG_IAC.Visible = 'off';
            app.FLAG_IAC.Value = true;
            gui_FLAG_IACValueChanged(app);
            % Update input items
            app.PP1.Visible = 'off'; app.PP2.Visible = 'off';
            app.PR3.Visible = 'on'; app.PR4.Visible = 'on';
            app.PP3.Visible = 'off'; app.PP4.Visible = 'off';
            app.PR5.Visible = 'off'; app.PP5.Visible = 'off';
            app.text_P1.Visible = 'on';
            app.text_RP2.Text = 'Pressure [bar]';
            % Visible flags
            app.FLAG_IAC.Visible = 'off';
            % Update Additional constraints panel
            app.AdditionalconstraintsPanel.Visible = 'on';
            app.AdditionalconstraintsPanel.Title = 'Additional constraints';
            app.text_RP.Visible ='off'; app.text_RP.Text = 'Products';
            app.text_R2.Visible = 'on'; app.text_P2.Visible = 'off';
            app.text_R2.Text = 'Reactants'; app.text_P2.Text = 'Products';
            app.text_RP3.Text = 'Shock velocity [m/s]'; 
            app.text_RP4.Visible = 'on'; app.text_RP4.Text = 'Mach number [-]';
            app.text_RP5.Visible = 'off'; 
            % Set default input values
            app.PR1.Value = '300';
            app.PR2.Value = '1';
            app.PR4.Value = '2';
            gui_compute_mach_or_velocity(app, 'Mach');
            % Set visible shock/detonation items
            gui_visible_shocks(app, true);
            % Set invisible oblique shock/detontions items
            gui_visible_oblique(app, false);
            % Set invisible rocket items
            gui_visible_rocket(app, false);
        case 'SHOCK_R' % * SHOCK_R: CALCULATE PLANAR POST-REFLECTED SHOCK STATE
            % Visible flags
            app.FLAG_IAC.Visible = 'off';
            app.FLAG_IAC.Value = true;
            gui_FLAG_IACValueChanged(app);
            % Update input items
            app.PP1.Visible = 'off'; app.PP2.Visible = 'off';
            app.PR3.Visible = 'on'; app.PR4.Visible = 'on';
            app.PP3.Visible = 'off'; app.PP4.Visible = 'off';
            app.PR5.Visible = 'off'; app.PP5.Visible = 'off';
            app.text_P1.Visible = 'on';
            % Visible flags
            app.FLAG_IAC.Visible = 'off';
            % Update Additional constraints panel
            app.AdditionalconstraintsPanel.Visible = 'on';
            app.AdditionalconstraintsPanel.Title = 'Additional constraints';
            app.text_RP.Visible ='off'; app.text_RP.Text = 'Products';
            app.text_R2.Visible = 'on'; app.text_P2.Visible = 'off';
            app.text_R2.Text = 'Reactants'; app.text_P2.Text = 'Products';
            app.text_RP3.Text = 'Shock velocity [m/s]'; 
            app.text_RP4.Visible = 'on'; app.text_RP4.Text = 'Mach number [-]';
            app.text_RP5.Visible = 'off';
            app.text_RP2.Text = 'Pressure [bar]';
            % Set default input values
            app.PR1.Value = '300';
            app.PR2.Value = '1';
            app.PR4.Value = '2';
            gui_compute_mach_or_velocity(app, 'Mach');
            % Set visible shocks/detonation items
            gui_visible_shocks(app, true);
            % Set invisible oblique shock/detontions items
            gui_visible_oblique(app, false);
            % Set invisible rocket items
            gui_visible_rocket(app, false);
        case {'SHOCK_OBLIQUE'}
            % Visible flags
            app.FLAG_IAC.Visible = 'off';
            app.FLAG_IAC.Value = true;
            gui_FLAG_IACValueChanged(app);
            % Update input items
            app.PP1.Visible = 'off'; app.PP2.Visible = 'off';
            app.PR3.Visible = 'on'; app.PR4.Visible = 'on';
            app.PP3.Visible = 'off'; app.PP4.Visible = 'off';
            app.PR5.Visible = 'on'; app.PP5.Visible = 'on';
            app.text_P1.Visible = 'on';
            app.text_RP2.Text = 'Pressure [bar]';
            % Visible flags
            app.FLAG_IAC.Visible = 'off';
            % Update Additional constraints panel
            app.AdditionalconstraintsPanel.Visible = 'on';
            app.AdditionalconstraintsPanel.Title = 'Additional constraints';
            app.text_RP.Visible ='on'; app.text_RP.Text = 'Reactants';
            app.text_R2.Visible = 'off'; app.text_P2.Visible = 'off';
            app.text_R2.Text = 'Reactants'; app.text_P2.Text = 'Products';
            app.text_RP3.Text = 'Shock velocity [m/s]'; 
            app.text_RP4.Visible = 'on'; app.text_RP4.Text = 'Mach number [-]';
            app.text_RP5.Visible = 'on';
            app.text_RP5.Text = 'Wave/Deflection angle [deg]';
            % Set default input values
            app.PR1.Value = '300';
            app.PR2.Value = '1';
            app.PR4.Value = '5';
            gui_compute_mach_or_velocity(app, 'Mach');
            app.PR5.Value = '40';
            app.PP5.Value = '';
            % Set visible shock/detonation items
            gui_visible_shocks(app, true);
            % Set visible oblique shock/detontions items
            gui_visible_oblique(app, true);
            % Set invisible rocket items
            gui_visible_rocket(app, false);
        case 'SHOCK_POLAR'
            % Visible flags
            app.FLAG_IAC.Visible = 'off';
            app.FLAG_IAC.Value = true;
            gui_FLAG_IACValueChanged(app);
            % Update input items
            app.PP1.Visible = 'off'; app.PP2.Visible = 'off';
            app.PR3.Visible = 'on'; app.PR4.Visible = 'on';
            app.PP3.Visible = 'off'; app.PP4.Visible = 'off';
            app.PR5.Visible = 'off'; app.PP5.Visible = 'off';
            app.text_P1.Visible = 'on';
            app.text_RP2.Text = 'Pressure [bar]';
            % Visible flags
            app.FLAG_IAC.Visible = 'off';
            % Update Additional constraints panel
            app.AdditionalconstraintsPanel.Visible = 'on';
            app.AdditionalconstraintsPanel.Title = 'Additional constraints';
            app.text_RP.Visible ='off'; app.text_RP.Text = 'Products';
            app.text_R2.Visible = 'on'; app.text_P2.Visible = 'off';
            app.text_R2.Text = 'Reactants'; app.text_P2.Text = 'Products';
            app.text_RP3.Text = 'Shock velocity [m/s]'; 
            app.text_RP4.Visible = 'on'; app.text_RP4.Text = 'Mach number [-]';
            app.text_RP5.Visible = 'off'; 
            % Set default input values
            app.PR1.Value = '300';
            app.PR2.Value = '1';
            app.PR4.Value = '2';
            gui_compute_mach_or_velocity(app, 'Mach');
            % Set visible shock/detonation items
            gui_visible_shocks(app, true);
            % Set invisible oblique shock/detontions items
            gui_visible_oblique(app, false);
            % Set invisible rocket items
            gui_visible_rocket(app, false);
        case 'SHOCK_POLAR_R'
            % Visible flags
            app.FLAG_IAC.Visible = 'off';
            app.FLAG_IAC.Value = true;
            gui_FLAG_IACValueChanged(app);
            % Update input items
            app.PP1.Visible = 'off'; app.PP2.Visible = 'off';
            app.PR3.Visible = 'on'; app.PR4.Visible = 'on';
            app.PP3.Visible = 'off'; app.PP4.Visible = 'off';
            app.PR5.Visible = 'on'; app.PP5.Visible = 'on';
            app.text_P1.Visible = 'on';
            app.text_RP2.Text = 'Pressure [bar]';
            % Visible flags
            app.FLAG_IAC.Visible = 'off';
            % Update Additional constraints panel
            app.AdditionalconstraintsPanel.Visible = 'on';
            app.AdditionalconstraintsPanel.Title = 'Additional constraints';
            app.text_RP.Visible ='on'; app.text_RP.Text = 'Reactants';
            app.text_R2.Visible = 'off'; app.text_P2.Visible = 'off';
            app.text_R2.Text = 'Reactants'; app.text_P2.Text = 'Products';
            app.text_RP3.Text = 'Shock velocity [m/s]'; 
            app.text_RP4.Visible = 'on'; app.text_RP4.Text = 'Mach number [-]';
            app.text_RP5.Visible = 'on';
            app.text_RP5.Text = 'Wave/Deflection angle [deg]';
            % Set default input values
            app.PR1.Value = '226.65';
            app.PR2.Value = '0.0117';
            app.PR4.Value = '20';
            gui_compute_mach_or_velocity(app, 'Mach');
            app.PR5.Value = '';
            app.PP5.Value = '35';
            % Set visible shock/detonation items
            gui_visible_shocks(app, true);
            % Set invisible oblique shock/detontions items
            gui_visible_oblique(app, false);
            % Set invisible rocket items
            gui_visible_rocket(app, false);
        case {'DET', 'DET_R'} % * DET: CALCULATE CHAPMAN-JOUGUET STATE
            % Visible flags
            app.FLAG_IAC.Visible = 'off';
            app.FLAG_IAC.Value = true;
            gui_FLAG_IACValueChanged(app);
            % Update input items
            app.PP1.Visible = 'off'; app.PP2.Visible = 'off';
            app.PR3.Visible = 'off'; app.PR4.Visible = 'off';
            app.PP3.Visible = 'off'; app.PP4.Visible = 'off';
            app.PR5.Visible = 'off'; app.PP5.Visible = 'off';
            app.text_P1.Visible = 'on';
            app.text_RP2.Text = 'Pressure [bar]';
            % Visible flags
            app.FLAG_IAC.Visible = 'off';
            % Update Additional constraints panel
            app.AdditionalconstraintsPanel.Visible = 'off';
            app.AdditionalconstraintsPanel.Title = 'Additional constraints';
            % Set default input values
            app.PR1.Value = '300';
            app.PR2.Value = '1';
            % Set visible shock/detonation items
            gui_visible_shocks(app, true);
            % Set invisible rocket items
            gui_visible_rocket(app, false);
        case {'DET_OVERDRIVEN', 'DET_OVERDRIVEN_R','DET_UNDERDRIVEN', 'DET_UNDERDRIVEN_R'} % * DET_OVERDRIVEN and DET_UNDERDRIVEN
            % Visible flags
            app.FLAG_IAC.Visible = 'off';
            app.FLAG_IAC.Value = true;
            gui_FLAG_IACValueChanged(app);
            app.PR5.Visible = 'off'; app.PP5.Visible = 'off';
            % Update input items
            app.PP1.Visible = 'off'; app.PP2.Visible = 'off';
            app.PR3.Visible = 'on'; app.PR4.Visible = 'off';
            app.PP3.Visible = 'off'; app.PP4.Visible = 'off';
            app.text_P1.Visible = 'on';
            app.text_RP2.Text = 'Pressure [bar]';
            % Visible flags
            app.FLAG_IAC.Visible = 'off';
            % Update Additional constraints panel
            app.AdditionalconstraintsPanel.Visible = 'on';
            app.AdditionalconstraintsPanel.Title = 'Additional constraints';
            app.text_RP.Visible ='off'; app.text_RP.Text = 'Products';
            app.text_R2.Visible = 'on'; app.text_P2.Visible = 'off';
            app.text_R2.Text = 'Reactants'; app.text_P2.Text = 'Products';
            if contains(app.ProblemType.Value, 'OVER')
                app.text_RP3.Text = 'Overdriven parameter [-]';
            else
                app.text_RP3.Text = 'Underdriven parameter [-]';
            end
            app.text_RP4.Visible = 'off'; 
            app.text_RP5.Visible = 'off'; 
            % Set default input values
            app.PR1.Value = '300';
            app.PR2.Value = '1';
            app.PR3.Value = '2';
            % Set visible shock/detonation items
            gui_visible_shocks(app, true);
            % Set invisible oblique shock/detontions items
            gui_visible_oblique(app, false);
            % Set invisible rocket items
            gui_visible_rocket(app, false);
        case {'DET_OBLIQUE'}
            % Visible flags
            app.FLAG_IAC.Visible = 'off';
            app.FLAG_IAC.Value = true;
            gui_FLAG_IACValueChanged(app);
            % Update input items
            app.PP1.Visible = 'off'; app.PP2.Visible = 'off';
            app.PR3.Visible = 'on'; app.PR4.Visible = 'on';
            app.PP3.Visible = 'off'; app.PP4.Visible = 'on';
            app.PR5.Visible = 'off'; app.PP5.Visible = 'off';
            app.text_P1.Visible = 'on';
            app.text_RP2.Text = 'Pressure [bar]';
            % Visible flags
            app.FLAG_IAC.Visible = 'off';
            % Update Additional constraints panel
            app.AdditionalconstraintsPanel.Visible = 'on';
            app.AdditionalconstraintsPanel.Title = 'Additional constraints';
            app.text_RP.Visible ='on'; app.text_RP.Text = 'Reactants';
            app.text_R2.Visible = 'off'; app.text_P2.Visible = 'off';
            app.text_R2.Text = 'Reactants'; app.text_P2.Text = 'Products';
            app.text_RP3.Text = 'Driven parameter [-]';
            app.text_RP4.Text = 'Wave/Deflection angle [deg]';
            app.text_RP4.Visible = 'on';
            app.text_RP5.Visible = 'off';
            % Set default input values
            app.PR1.Value = '300';
            app.PR2.Value = '1';
            app.PR3.Value = '2';
            app.PR4.Value = '60';
            app.PP4.Value = '';
            % Set visible shock/detonation items
            gui_visible_shocks(app, true);
            % Set visible oblique shock/detontions items
            gui_visible_oblique(app, true);
            % Set invisible rocket items
            gui_visible_rocket(app, false);
        case {'DET_POLAR'}
            % Visible flags
            app.FLAG_IAC.Visible = 'off';
            app.FLAG_IAC.Value = true;
            gui_FLAG_IACValueChanged(app);
            app.PR5.Visible = 'off'; app.PP5.Visible = 'off';
            % Update input items
            app.PP1.Visible = 'off'; app.PP2.Visible = 'off';
            app.PR3.Visible = 'on'; app.PR4.Visible = 'off';
            app.PP3.Visible = 'off'; app.PP4.Visible = 'off';
            app.text_P1.Visible = 'on';
            app.text_RP2.Text = 'Pressure [bar]';
            % Visible flags
            app.FLAG_IAC.Visible = 'off';
            % Update Additional constraints panel
            app.AdditionalconstraintsPanel.Visible = 'on';
            app.AdditionalconstraintsPanel.Title = 'Additional constraints';
            app.text_RP.Visible ='off'; app.text_RP.Text = 'Products';
            app.text_R2.Visible = 'on'; app.text_P2.Visible = 'off';
            app.text_R2.Text = 'Reactants'; app.text_P2.Text = 'Products';
            app.text_RP3.Text = 'Overdriven parameter [-]';
            app.text_RP4.Visible = 'off'; 
            app.text_RP5.Visible = 'off'; 
            % Set default input values
            app.PR1.Value = '300';
            app.PR2.Value = '1';
            app.PR3.Value = '2';
            % Set visible shock/detonation items
            gui_visible_shocks(app, true);
            % Set invisible oblique shock/detontions items
            gui_visible_oblique(app, false);
            % Set invisible rocket items
            gui_visible_rocket(app, false);
        case {'ROCKET'}
            % Visible flags
            app.FLAG_IAC.Visible = 'on';
            app.FLAG_IAC.Value = true;
            gui_FLAG_IACValueChanged(app);
            % Update input items
            app.PP1.Visible = 'off'; app.PP2.Visible = 'off';
            app.PR3.Visible = 'on'; app.PR4.Visible = 'off';
            app.PP3.Visible = 'on'; app.PP4.Visible = 'off';
            app.PR5.Visible = 'off'; app.PP5.Visible = 'off';
            app.text_P1.Visible = 'off';
            app.text_RP2.Text = 'Pressure [bar]';
            % Update Additional constraints panel
            app.AdditionalconstraintsPanel.Visible = 'on';
            app.AdditionalconstraintsPanel.Title = 'Optional parameters';
            app.text_RP3.Text = 'Area ratios A/A_throat';
            app.text_RP.Visible = 'off'; app.text_RP.Text = 'Products';
            app.text_R2.Visible = 'on'; app.text_P2.Visible = 'on';
            app.text_R2.Text = 'Subsonic'; app.text_P2.Text = 'Supersonic';
            app.text_RP4.Visible = 'off';
            app.text_RP5.Visible = 'off'; 
            % Set default input values
            app.PR1.Value = '300';
            app.PR2.Value = '1';
            app.PR3.Value = '';
            app.PP3.Value = '';
            % Set visible shock/detonation items
            gui_visible_shocks(app, true);
            % Set invisible oblique shock/detontions items
            gui_visible_oblique(app, false);
            % Set visible rocket items
            gui_visible_rocket(app, true);
         case {'SHOCK_OBLIQUE_R', 'DET_OBLIQUE_R', 'DET_POLAR_R'}
            uialert(app.UIFigure, 'Problem not included yet. Sorry for the inconvinience.', 'Error')
    end

    % Update temperature 
    if isempty(app.UITable_R.Data)
        return
    end
    
    app.UITable_R.Data(:, 5) = repmat({gui_get_prop(app.PR1.Value, 'first')}, size(app.UITable_R.Data(:, 5)));
end

% SUB-PASS FUNCTIONS
function gui_visible_shocks(app, value)
    if value
        value = 'on';
    else
        value = 'off';
    end
    app.text_u.Visible = value; app.text_u_1.Visible = value; app.text_u_2.Visible = value;
    app.text_M.Visible = value; app.text_M_1.Visible = value; app.text_M_2.Visible = value;
end

function gui_visible_oblique(app, value)
    if value
        value = 'on';
    else
        value = 'off';
    end
    app.text_beta_min.Visible = value; app.text_beta_min_2.Visible = value;
    app.text_beta.Visible = value; app.text_beta_2.Visible = value;
    app.text_theta.Visible = value; app.text_theta_2.Visible = value;
end

function gui_visible_rocket(app, value)
    if value
        value = 'on';
        delta_pos_x = 77;
        % Move panel properties
        app.Panel_parameters.Position(1) = app.default.Panel_parameters.Position(1) - delta_pos_x;
        app.Panel_parameters.Position(3) = app.default.Panel_parameters.Position(3) + 3 * delta_pos_x;
        % Change GUI to default rocket model
        app.public_FLAG_IACValueChanged();
    else
        value = 'off';
        app.Panel_parameters.Position = app.default.Panel_parameters.Position;
        app.text_Products.Text = 'Products';
        app.Panel_extra_5.Visible = 'off';
    end
    % Set visible properties of Rockets
    app.text_Aratio.Visible = value; app.text_Aratio_2.Visible = value;
    app.text_Cstar.Visible = value;
    app.text_Ivac.Visible = value;
    app.text_Isp.Visible = value;
    app.Panel_extra_3.Visible = value;
    app.Panel_extra_4.Visible = value;
end