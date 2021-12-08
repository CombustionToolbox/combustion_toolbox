function gui_ProblemTypeValueChanged(obj)
    % Clear GUI results tab (except UITree) and update GUI items for the problem selected

    % 1. Clear results tab (except UITree)
    gui_clear_results(obj);

    % 2. Update GUI items depending of the problem selected    
    switch obj.ProblemType.Value
        case 'TP' % * TP: Equilibrium composition at defined T and p
            % Update input items
            obj.PP1.Visible = 'on'; obj.PP2.Visible = 'on';
            % Update Additional constraints panel
            obj.AdditionalconstraintsPanel.Visible = 'off';
            % Set default input values
            obj.PR1.Value = '300'; 
            obj.PR2.Value = '1';
            obj.PP1.Value = '2500';
            obj.PP2.Value = obj.PR2.Value;
            % Set visible item volume ratio
            obj.text_vP_vR.Visible = 'off'; 
            % Set invisible shocks/detonation items
            gui_visible_shocks(obj, false);
        case 'HP' % * HP: Adiabatic T and composition at constant p
            % Update input items
            obj.PP1.Visible = 'off'; obj.PP2.Visible = 'on';
            obj.PR3.Visible = 'off'; obj.PR4.Visible = 'off';
            obj.PP3.Visible = 'off'; obj.PP4.Visible = 'off'; 
            % Update Additional constraints panel
            obj.AdditionalconstraintsPanel.Visible = 'on';
            obj.text_RP3.Text = 'Constant Enthalpy: hP = hR';
            obj.text_RP.Visible = 'on'; obj.text_R2.Visible = 'off'; obj.text_P2.Visible = 'off';
            obj.text_RP4.Visible = 'off';
            % Set default input values
            obj.PR1.Value = '300';
            obj.PR2.Value = '1';
            obj.PP2.Value = obj.PR2.Value;
            % Set visible item volume ratio
            obj.text_vP_vR.Visible = 'off'; 
            % Set invisible shocks/detonation items
            gui_visible_shocks(obj, false);
        case 'SP' % * SP: Isentropic (i.e., adiabatic) compression/expansion to a specified p
            % Update input items
            obj.PP1.Visible = 'off'; obj.PP2.Visible = 'on';
            obj.PR3.Visible = 'off'; obj.PR4.Visible = 'off';
            obj.PP3.Visible = 'off'; obj.PP4.Visible = 'off';
            % Update Additional constraints panel
            obj.AdditionalconstraintsPanel.Visible = 'on';
            obj.text_RP.Visible = 'on'; obj.text_R2.Visible = 'off'; obj.text_P2.Visible = 'off';
            obj.text_RP3.Text = 'Constant Entropy: SP = SR';
            obj.text_RP4.Visible = 'off'; 
            % Set default input values
            obj.PR1.Value = '300';
            obj.PR2.Value = '20';
            obj.PP2.Value = obj.PR2.Value;
            % Set visible item volume ratio
            obj.text_vP_vR.Visible = 'off'; 
            % Set invisible shocks/detonation items
            gui_visible_shocks(obj, false);
        case 'TV' % * TV: Equilibrium composition at defined T and constant v
            % Update input items
            obj.PP1.Visible = 'on'; obj.PP2.Visible = 'off'; 
            obj.PR3.Visible = 'off'; obj.PR4.Visible = 'off';
            obj.PP3.Visible = 'off'; obj.PP4.Visible = 'off'; 
            % Update Additional constraints panel
            obj.AdditionalconstraintsPanel.Visible = 'on';
            obj.text_RP.Visible = 'on'; obj.text_R2.Visible = 'off'; obj.text_P2.Visible = 'off';
            obj.text_RP3.Text = 'Constant Volume: vP = vR';
            obj.text_RP4.Visible = 'off'; 
            % Set default input values
            obj.PR1.Value = '300';
            obj.PR2.Value = '1';
            obj.PP1.Value = '2500';
            obj.PP2.Value = obj.PR2.Value;
            % Set visible item volume ratio
            obj.text_vP_vR.Visible = 'off'; 
            % Set invisible shocks/detonation items
            gui_visible_shocks(obj, false);
        case 'EV' % * EV: Equilibrium composition at Adiabatic T and constant v
            % Update input items
            obj.PP1.Visible = 'off'; obj.PP2.Visible = 'off'; 
            obj.PR3.Visible = 'off'; obj.PR4.Visible = 'off';
            obj.PP3.Visible = 'off'; obj.PP4.Visible = 'off'; 
            % Update Additional constraints panel
            obj.AdditionalconstraintsPanel.Visible = 'on';
            obj.text_RP.Visible = 'on'; obj.text_R2.Visible = 'off'; obj.text_P2.Visible = 'off';
            obj.text_RP3.Text = 'Constant Internal energy: eP = eR';
            obj.text_RP4.Visible = 'on'; obj.text_RP4.Text = 'Constant Volume: vP = vR';
            % Set default input values
            obj.PR1.Value = '1000';
            obj.PR2.Value = '1';
            obj.PP2.Value = obj.PR2.Value;
            % Set visible item volume ratio
            obj.text_vP_vR.Visible = 'off'; 
            % Set invisible shocks/detonation items
            gui_visible_shocks(obj, false);
        case 'SV' % * SV: Isentropic (i.e., fast adiabatic) compression/expansion to a specified v
            % Update input items
            obj.PP1.Visible = 'off'; obj.PP2.Visible = 'off'; 
            obj.PR3.Visible = 'off'; obj.PR4.Visible = 'off';
            obj.PP3.Visible = 'off'; obj.PP4.Visible = 'on'; 
            % Update Additional constraints panel
            obj.AdditionalconstraintsPanel.Visible = 'on';
            obj.text_RP.Visible = 'on'; obj.text_R2.Visible = 'off'; obj.text_P2.Visible = 'off';
            obj.text_RP3.Text = 'Constant Entropy: SP = SR';
            obj.text_RP4.Visible = 'on'; obj.text_RP4.Text = 'Volume Products/Reactants';
            % Set default input values
            obj.PR1.Value = '300';
            obj.PR2.Value = '1';
            obj.PP2.Value = obj.PR2.Value;
            obj.PP4.Value = '0.5';
            % Set visible item volume ratio
            obj.text_vP_vR.Visible = 'on'; 
            % Set invisible shocks/detonation items
            gui_visible_shocks(obj, false);
        case 'SHOCK_I' % * SHOCK_I: CALCULATE PLANAR INCIDENT SHOCK WAVE
            % Update input items
            obj.PP1.Visible = 'off'; obj.PP2.Visible = 'off';
            obj.PR3.Visible = 'on'; obj.PR4.Visible = 'on';
            obj.PP3.Visible = 'off'; obj.PP4.Visible = 'off'; 
            % Update Additional constraints panel
            obj.AdditionalconstraintsPanel.Visible = 'on';
            obj.text_RP.Visible ='off'; obj.text_R2.Visible = 'on'; obj.text_P2.Visible = 'off';
            obj.text_RP3.Text = 'Shock velocity [m/s]'; 
            obj.text_RP4.Visible = 'on'; obj.text_RP4.Text = 'Mach number [-]'; 
            % Set default input values
            obj.PR1.Value = '300';
            obj.PR2.Value = '1';
            obj.PR4.Value = '2';
            gui_compute_mach_or_velocity(obj, 'Mach');
            % Set visible item volume ratio
            obj.text_vP_vR.Visible = 'off';
            % Set invisible shocks/detonation items
            gui_visible_shocks(obj, true);
        case 'SHOCK_R' % * SHOCK_R: CALCULATE PLANAR POST-REFLECTED SHOCK STATE
            % Update input items
            obj.PP1.Visible = 'off'; obj.PP2.Visible = 'off';
            obj.PR3.Visible = 'on'; obj.PR4.Visible = 'on';
            obj.PP3.Visible = 'off'; obj.PP4.Visible = 'off'; 
            % Update Additional constraints panel
            obj.AdditionalconstraintsPanel.Visible = 'on';
            obj.text_RP.Visible ='off'; obj.text_R2.Visible = 'on'; obj.text_P2.Visible = 'off';
            obj.text_RP3.Text = 'Shock velocity [m/s]'; 
            obj.text_RP4.Visible = 'on'; obj.text_RP4.Text = 'Mach number [-]';
            % Set default input values
            obj.PR1.Value = '300';
            obj.PR2.Value = '1';
            obj.PR4.Value = '2';
            gui_compute_mach_or_velocity(obj, 'Mach');
            % Set visible item volume ratio
            obj.text_vP_vR.Visible = 'off'; 
            % Set invisible shocks/detonation items
            gui_visible_shocks(obj, true);
        case 'DET' % * DET: CALCULATE CHAPMAN-JOUGET STATE (CJ UPPER STATE)
            % Update input items
            obj.PP1.Visible = 'off'; obj.PP2.Visible = 'off';
            obj.PR3.Visible = 'off'; obj.PR4.Visible = 'off';
            obj.PP3.Visible = 'off'; obj.PP4.Visible = 'off'; 
            % Update Additional constraints panel
            obj.AdditionalconstraintsPanel.Visible = 'off';
            % Set default input values
            obj.PR1.Value = '300';
            obj.PR2.Value = '1';
            % Set visible item volume ratio
            obj.text_vP_vR.Visible = 'off'; 
            % Set invisible shocks/detonation items
            gui_visible_shocks(obj, true);
        case 'DET_OVERDRIVEN' % * DET_OVERDRIVEN: CALCULATE OVERDRIVE DETONATION
            % Update input items
            obj.PP1.Visible = 'off'; obj.PP2.Visible = 'off';
            obj.PR3.Visible = 'on'; obj.PR4.Visible = 'off';
            obj.PP3.Visible = 'off'; obj.PP4.Visible = 'off'; 
            % Update Additional constraints panel
            obj.AdditionalconstraintsPanel.Visible = 'on';
            obj.text_RP.Visible ='off'; obj.text_R2.Visible = 'on'; obj.text_P2.Visible = 'off';
            obj.text_RP3.Text = 'Overdriven parameter [-]'; 
            obj.text_RP4.Visible = 'off'; 
            % Set default input values
            obj.PR1.Value = '300';
            obj.PR2.Value = '1';
            obj.PR3.Value = '2';
            % Set visible item volume ratio
            obj.text_vP_vR.Visible = 'off'; 
            % Set invisible shocks/detonation items
            gui_visible_shocks(obj, true);
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