function gui_ProblemTypeValueChanged(obj)
    % Clear GUI results tab (except UITree) and update GUI items for the problem selected

    % 1. Clear results tab (except UITree)
    gui_clear_results(obj);

    % 2. Update GUI items depending of the problem selected    
    switch obj.ProblemType.Value
        case 'TP' % * TP: Equilibrium composition at defined T and p
            % Visible flags
            obj.FLAG_IAC.Visible = 'off';
            obj.FLAG_IAC.Value = true;
            gui_FLAG_IACValueChanged(obj);
            % Update input items
            obj.PP1.Visible = 'on'; obj.PP2.Visible = 'on';
            obj.text_P1.Visible = 'on';
            % Update Additional constraints panel
            obj.AdditionalconstraintsPanel.Visible = 'off';
            obj.AdditionalconstraintsPanel.Title = 'Additional constraints';
            % Set default input values
            obj.PR1.Value = '300'; 
            obj.PR2.Value = '1';
            obj.PP1.Value = '2500';
            obj.PP2.Value = obj.PR2.Value;
            % Set invisible shocks/detonation items
            gui_visible_shocks(obj, false);
            % Set invisible rocket items
            gui_visible_rocket(obj, false);
        case 'HP' % * HP: Adiabatic T and composition at constant p
            % Visible flags
            obj.FLAG_IAC.Visible = 'off';
            obj.FLAG_IAC.Value = true;
            gui_FLAG_IACValueChanged(obj);
            % Update input items
            obj.PP1.Visible = 'off'; obj.PP2.Visible = 'on';
            obj.PR3.Visible = 'off'; obj.PR4.Visible = 'off';
            obj.PP3.Visible = 'off'; obj.PP4.Visible = 'off'; 
            obj.text_P1.Visible = 'on';
            % Visible flags
            obj.FLAG_IAC.Visible = 'off';
            % Update Additional constraints panel
            obj.AdditionalconstraintsPanel.Visible = 'on';
            obj.AdditionalconstraintsPanel.Title = 'Additional constraints';
            obj.text_RP3.Text = 'Constant Enthalpy: hP = hR';
            obj.text_RP.Visible = 'on'; obj.text_R2.Visible = 'off'; obj.text_P2.Visible = 'off';
            obj.text_R2.Text = 'Reactants'; obj.text_P2.Text = 'Products';
            obj.text_RP4.Visible = 'off';
            % Set default input values
            obj.PR1.Value = '300';
            obj.PR2.Value = '1';
            obj.PP2.Value = obj.PR2.Value;
            % Set invisible shocks/detonation items
            gui_visible_shocks(obj, false);
            % Set invisible rocket items
            gui_visible_rocket(obj, false);
        case 'SP' % * SP: Isentropic (i.e., adiabatic) compression/expansion to a specified p
            % Visible flags
            obj.FLAG_IAC.Visible = 'off';
            obj.FLAG_IAC.Value = true;
            gui_FLAG_IACValueChanged(obj);
            % Update input items
            obj.PP1.Visible = 'off'; obj.PP2.Visible = 'on';
            obj.PR3.Visible = 'off'; obj.PR4.Visible = 'off';
            obj.PP3.Visible = 'off'; obj.PP4.Visible = 'off';
            obj.text_P1.Visible = 'on';
            % Visible flags
            obj.FLAG_IAC.Visible = 'off';
            % Update Additional constraints panel
            obj.AdditionalconstraintsPanel.Visible = 'on';
            obj.AdditionalconstraintsPanel.Title = 'Additional constraints';
            obj.text_RP.Visible = 'on'; obj.text_R2.Visible = 'off'; obj.text_P2.Visible = 'off';
            obj.text_R2.Text = 'Reactants'; obj.text_P2.Text = 'Products';
            obj.text_RP3.Text = 'Constant Entropy: SP = SR';
            obj.text_RP4.Visible = 'off'; 
            % Set default input values
            obj.PR1.Value = '300';
            obj.PR2.Value = '1';
            obj.PP2.Value = '100';
            % Set invisible shocks/detonation items
            gui_visible_shocks(obj, false);
            % Set invisible rocket items
            gui_visible_rocket(obj, false);
        case 'TV' % * TV: Equilibrium composition at defined T and constant v
            % Visible flags
            obj.FLAG_IAC.Visible = 'off';
            obj.FLAG_IAC.Value = true;
            gui_FLAG_IACValueChanged(obj);
            % Update input items
            obj.PP1.Visible = 'on'; obj.PP2.Visible = 'off'; 
            obj.PR3.Visible = 'off'; obj.PR4.Visible = 'off';
            obj.PP3.Visible = 'off'; obj.PP4.Visible = 'off'; 
            obj.text_P1.Visible = 'on';
            % Visible flags
            obj.FLAG_IAC.Visible = 'off';
            % Update Additional constraints panel
            obj.AdditionalconstraintsPanel.Visible = 'on';
            obj.AdditionalconstraintsPanel.Title = 'Additional constraints';
            obj.text_RP.Visible = 'on'; obj.text_R2.Visible = 'off'; obj.text_P2.Visible = 'off';
            obj.text_R2.Text = 'Reactants'; obj.text_P2.Text = 'Products';
            obj.text_RP3.Text = 'Constant Volume: vP = vR';
            obj.text_RP4.Visible = 'off'; 
            % Set default input values
            obj.PR1.Value = '300';
            obj.PR2.Value = '1';
            obj.PP1.Value = '2500';
            obj.PP2.Value = obj.PR2.Value;
            % Set invisible shocks/detonation items
            gui_visible_shocks(obj, false);
            % Set invisible rocket items
            gui_visible_rocket(obj, false);
        case 'EV' % * EV: Equilibrium composition at Adiabatic T and constant v
            % Visible flags
            obj.FLAG_IAC.Visible = 'off';
            obj.FLAG_IAC.Value = true;
            gui_FLAG_IACValueChanged(obj);
            % Update input items
            obj.PP1.Visible = 'off'; obj.PP2.Visible = 'off'; 
            obj.PR3.Visible = 'off'; obj.PR4.Visible = 'off';
            obj.PP3.Visible = 'off'; obj.PP4.Visible = 'off'; 
            obj.text_P1.Visible = 'on';
            % Visible flags
            obj.FLAG_IAC.Visible = 'off';
            % Update Additional constraints panel
            obj.AdditionalconstraintsPanel.Visible = 'on';
            obj.AdditionalconstraintsPanel.Title = 'Additional constraints';
            obj.text_RP.Visible = 'on'; obj.text_R2.Visible = 'off'; obj.text_P2.Visible = 'off';
            obj.text_R2.Text = 'Reactants'; obj.text_P2.Text = 'Products';
            obj.text_RP3.Text = 'Constant Internal energy: eP = eR';
            obj.text_RP4.Visible = 'on'; obj.text_RP4.Text = 'Constant Volume: vP = vR';
            % Set default input values
            obj.PR1.Value = '1000';
            obj.PR2.Value = '1';
            obj.PP2.Value = obj.PR2.Value;
            % Set invisible shocks/detonation items
            gui_visible_shocks(obj, false);
            % Set invisible rocket items
            gui_visible_rocket(obj, false);
        case 'SV' % * SV: Isentropic (i.e., fast adiabatic) compression/expansion to a specified v
            % Visible flags
            obj.FLAG_IAC.Visible = 'off';
            obj.FLAG_IAC.Value = true;
            gui_FLAG_IACValueChanged(obj);
            % Update input items
            obj.PP1.Visible = 'off'; obj.PP2.Visible = 'off'; 
            obj.PR3.Visible = 'off'; obj.PR4.Visible = 'off';
            obj.PP3.Visible = 'off'; obj.PP4.Visible = 'on'; 
            obj.text_P1.Visible = 'on';
            % Visible flags
            obj.FLAG_IAC.Visible = 'off';
            % Update Additional constraints panel
            obj.AdditionalconstraintsPanel.Visible = 'on';
            obj.AdditionalconstraintsPanel.Title = 'Additional constraints';
            obj.text_RP.Visible = 'on'; obj.text_R2.Visible = 'off'; obj.text_P2.Visible = 'off';
            obj.text_R2.Text = 'Reactants'; obj.text_P2.Text = 'Products';
            obj.text_RP3.Text = 'Constant Entropy: SP = SR';
            obj.text_RP4.Visible = 'on'; obj.text_RP4.Text = 'Volume Products/Reactants';
            % Set default input values
            obj.PR1.Value = '300';
            obj.PR2.Value = '1';
            obj.PP2.Value = obj.PR2.Value;
            % Set invisible shocks/detonation items
            gui_visible_shocks(obj, false);
            % Set invisible rocket items
            gui_visible_rocket(obj, false);
        case 'SHOCK_I' % * SHOCK_I: CALCULATE PLANAR INCIDENT SHOCK WAVE
            % Visible flags
            obj.FLAG_IAC.Visible = 'off';
            obj.FLAG_IAC.Value = true;
            gui_FLAG_IACValueChanged(obj);
            % Update input items
            obj.PP1.Visible = 'off'; obj.PP2.Visible = 'off';
            obj.PR3.Visible = 'on'; obj.PR4.Visible = 'on';
            obj.PP3.Visible = 'off'; obj.PP4.Visible = 'off'; 
            obj.text_P1.Visible = 'on';
            % Visible flags
            obj.FLAG_IAC.Visible = 'off';
            % Update Additional constraints panel
            obj.AdditionalconstraintsPanel.Visible = 'on';
            obj.AdditionalconstraintsPanel.Title = 'Additional constraints';
            obj.text_RP.Visible ='off'; obj.text_R2.Visible = 'on'; obj.text_P2.Visible = 'off';
            obj.text_R2.Text = 'Reactants'; obj.text_P2.Text = 'Products';
            obj.text_RP3.Text = 'Shock velocity [m/s]'; 
            obj.text_RP4.Visible = 'on'; obj.text_RP4.Text = 'Mach number [-]'; 
            % Set default input values
            obj.PR1.Value = '300';
            obj.PR2.Value = '1';
            obj.PR4.Value = '2';
            gui_compute_mach_or_velocity(obj, 'Mach');
            % Set visible shocks/detonation items
            gui_visible_shocks(obj, true);
        case 'SHOCK_R' % * SHOCK_R: CALCULATE PLANAR POST-REFLECTED SHOCK STATE
            % Visible flags
            obj.FLAG_IAC.Visible = 'off';
            obj.FLAG_IAC.Value = true;
            gui_FLAG_IACValueChanged(obj);
            % Update input items
            obj.PP1.Visible = 'off'; obj.PP2.Visible = 'off';
            obj.PR3.Visible = 'on'; obj.PR4.Visible = 'on';
            obj.PP3.Visible = 'off'; obj.PP4.Visible = 'off'; 
            obj.text_P1.Visible = 'on';
            % Visible flags
            obj.FLAG_IAC.Visible = 'off';
            % Update Additional constraints panel
            obj.AdditionalconstraintsPanel.Visible = 'on';
            obj.AdditionalconstraintsPanel.Title = 'Additional constraints';
            obj.text_RP.Visible ='off'; obj.text_R2.Visible = 'on'; obj.text_P2.Visible = 'off';
            obj.text_R2.Text = 'Reactants'; obj.text_P2.Text = 'Products';
            obj.text_RP3.Text = 'Shock velocity [m/s]'; 
            obj.text_RP4.Visible = 'on'; obj.text_RP4.Text = 'Mach number [-]';
            % Set default input values
            obj.PR1.Value = '300';
            obj.PR2.Value = '1';
            obj.PR4.Value = '2';
            gui_compute_mach_or_velocity(obj, 'Mach');
            % Set visible shocks/detonation items
            gui_visible_shocks(obj, true);
            % Set invisible rocket items
            gui_visible_rocket(obj, false);
        case {'DET', 'DET_R'} % * DET: CALCULATE CHAPMAN-JOUGUET STATE
            % Visible flags
            obj.FLAG_IAC.Visible = 'off';
            obj.FLAG_IAC.Value = true;
            gui_FLAG_IACValueChanged(obj);
            % Update input items
            obj.PP1.Visible = 'off'; obj.PP2.Visible = 'off';
            obj.PR3.Visible = 'off'; obj.PR4.Visible = 'off';
            obj.PP3.Visible = 'off'; obj.PP4.Visible = 'off'; 
            obj.text_P1.Visible = 'on';
            % Visible flags
            obj.FLAG_IAC.Visible = 'off';
            % Update Additional constraints panel
            obj.AdditionalconstraintsPanel.Visible = 'off';
            obj.AdditionalconstraintsPanel.Title = 'Additional constraints';
            % Set default input values
            obj.PR1.Value = '300';
            obj.PR2.Value = '1';
            % Set visible shocks/detonation items
            gui_visible_shocks(obj, true);
            % Set invisible rocket items
            gui_visible_rocket(obj, false);
        case {'DET_OVERDRIVEN', 'DET_OVERDRIVEN_R','DET_UNDERDRIVEN', 'DET_UNDERDRIVEN_R'} % * DET_OVERDRIVEN and DET_UNDERDRIVEN
            % Visible flags
            obj.FLAG_IAC.Visible = 'off';
            obj.FLAG_IAC.Value = true;
            gui_FLAG_IACValueChanged(obj);
            % Update input items
            obj.PP1.Visible = 'off'; obj.PP2.Visible = 'off';
            obj.PR3.Visible = 'on'; obj.PR4.Visible = 'off';
            obj.PP3.Visible = 'off'; obj.PP4.Visible = 'off';
            obj.text_P1.Visible = 'on';
            % Visible flags
            obj.FLAG_IAC.Visible = 'off';
            % Update Additional constraints panel
            obj.AdditionalconstraintsPanel.Visible = 'on';
            obj.AdditionalconstraintsPanel.Title = 'Additional constraints';
            obj.text_RP.Visible ='off'; obj.text_R2.Visible = 'on'; obj.text_P2.Visible = 'off';
            obj.text_R2.Text = 'Reactants'; obj.text_P2.Text = 'Products';
            if contains(obj.ProblemType.Value, 'OVER')
                obj.text_RP3.Text = 'Overdriven parameter [-]';
            else
                obj.text_RP3.Text = 'Underdriven parameter [-]';
            end
            obj.text_RP4.Visible = 'off'; 
            % Set default input values
            obj.PR1.Value = '300';
            obj.PR2.Value = '1';
            obj.PR3.Value = '2';
            % Set visible shocks/detonation items
            gui_visible_shocks(obj, true);
            % Set invisible rocket items
            gui_visible_rocket(obj, false);
        case {'ROCKET'} % * ROCKET: ROCKET PROPELLANT PERFORMANCE
            % Visible flags
            obj.FLAG_IAC.Visible = 'on';
            obj.FLAG_IAC.Value = true;
            gui_FLAG_IACValueChanged(obj);
            % Update input items
            obj.PP1.Visible = 'off'; obj.PP2.Visible = 'off';
            obj.PR3.Visible = 'on'; obj.PR4.Visible = 'off';
            obj.PP3.Visible = 'on'; obj.PP4.Visible = 'off'; 
            obj.text_P1.Visible = 'off';
            % Update Additional constraints panel
            obj.AdditionalconstraintsPanel.Visible = 'on';
            obj.AdditionalconstraintsPanel.Title = 'Optional parameters';
            obj.text_RP3.Text = 'Area ratios A/A_throat';
            obj.text_RP.Visible = 'off'; obj.text_R2.Visible = 'on'; obj.text_P2.Visible = 'on';
            obj.text_R2.Text = 'Subsonic'; obj.text_P2.Text = 'Supersonic';
            obj.text_RP4.Visible = 'off';
            % Set default input values
            obj.PR1.Value = '300';
            obj.PR2.Value = '1';
            obj.PR3.Value = '';
            obj.PP3.Value = '';
            % Set visible shocks/detonation items
            gui_visible_shocks(obj, true);
            % Set visible rocket items
            gui_visible_rocket(obj, true);
    end
end

% SUB-PASS FUNCTIONS
function gui_visible_shocks(obj, value)
    if value
        value = 'on';
    else
        value = 'off';
    end
    obj.text_u.Visible = value; obj.text_uR.Visible = value; obj.text_uP.Visible = value;
    obj.text_MR.Visible = value; obj.text_MP.Visible = value; obj.text_M.Visible = value;
end

function gui_visible_rocket(obj, value)
    if value
        value = 'on';
        delta_pos_x = 77;
        % Move panel properties
        obj.Panel_parameters.Position(1) = obj.default.Panel_parameters.Position(1) - delta_pos_x;
        obj.Panel_parameters.Position(3) = obj.default.Panel_parameters.Position(3) + 3 * delta_pos_x;
        % Change GUI to default rocket model
        obj.public_FLAG_IACValueChanged();
    else
        value = 'off';
        obj.Panel_parameters.Position = obj.default.Panel_parameters.Position;
        obj.text_Products_2.Text = 'Products';
        obj.Panel_extra_5.Visible = 'off';
    end
    % Set visible properties of Rockets
    obj.text_Aratio.Visible = value; obj.text_Aratio_2.Visible = value;
    obj.text_Cstar.Visible = value;
    obj.text_Ivac.Visible = value;
    obj.text_Isp.Visible = value;
    obj.Panel_extra_3.Visible = value;
    obj.Panel_extra_4.Visible = value;
end