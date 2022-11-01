function gui_generate_panel_base(obj)
    % Routine that generates all the components to read the thermodynamic
    % properties at different stages of the rocket
    
    % Create hPhRkJkgEditFieldLabel
    obj.hPhRkJkgEditFieldLabel = uilabel(obj.ParametersTab);
    obj.hPhRkJkgEditFieldLabel.HorizontalAlignment = 'right';
    obj.hPhRkJkgEditFieldLabel.Position = [103 46 84 22];
    obj.hPhRkJkgEditFieldLabel.Text = 'hP - hR [kJ/kg]';

    % Create text_q
    obj.text_q = uieditfield(obj.ParametersTab, 'numeric');
    obj.text_q.Editable = 'off';
    obj.text_q.Position = [202 46 91 19];

    % Create VolumeProductsReactantsLabel
    obj.VolumeProductsReactantsLabel = uilabel(obj.ParametersTab);
    obj.VolumeProductsReactantsLabel.HorizontalAlignment = 'right';
    obj.VolumeProductsReactantsLabel.Position = [33 21 154 22];
    obj.VolumeProductsReactantsLabel.Text = 'Volume Products/Reactants';

    % Create text_vP_vR
    obj.text_vP_vR = uieditfield(obj.ParametersTab, 'numeric');
    obj.text_vP_vR.Editable = 'off';
    obj.text_vP_vR.Position = [202 21 91 19];

    % Create EpsilonmolesLabel
    obj.EpsilonmolesLabel = uilabel(obj.ParametersTab);
    obj.EpsilonmolesLabel.HorizontalAlignment = 'right';
    obj.EpsilonmolesLabel.Position = [317 46 88 22];
    obj.EpsilonmolesLabel.Text = 'Epsilon (moles)';

    % Create text_error_moles
    obj.text_error_moles = uieditfield(obj.ParametersTab, 'numeric');
    obj.text_error_moles.ValueDisplayFormat = '%11.4e';
    obj.text_error_moles.ValueChangedFcn = createCallbackFcn(obj, @text_error_molesValueChanged, true);
    obj.text_error_moles.Editable = 'off';
    obj.text_error_moles.Position = [420 46 91 19];

    % Create EpsilonmolesLabel_2
    obj.EpsilonmolesLabel_2 = uilabel(obj.ParametersTab);
    obj.EpsilonmolesLabel_2.HorizontalAlignment = 'right';
    obj.EpsilonmolesLabel_2.Position = [309 21 96 22];
    obj.EpsilonmolesLabel_2.Text = 'Epsilon (method)';

    % Create text_error_problem
    obj.text_error_problem = uieditfield(obj.ParametersTab, 'numeric');
    obj.text_error_problem.ValueDisplayFormat = '%11.4e';
    obj.text_error_problem.ValueChangedFcn = createCallbackFcn(obj, @text_error_problemValueChanged, true);
    obj.text_error_problem.Editable = 'off';
    obj.text_error_problem.Position = [420 21 91 19];

    % Create Panel_parameters
    obj.Panel_parameters = uipanel(obj.ParametersTab);
    obj.Panel_parameters.AutoResizeChildren = 'off';
    obj.Panel_parameters.BorderType = 'none';
    obj.Panel_parameters.BackgroundColor = [0.9098 0.9098 0.8902];
    obj.Panel_parameters.Position = [11 86 549 357];

    % Create text_s
    obj.text_s = uilabel(obj.Panel_parameters);
    obj.text_s.HorizontalAlignment = 'center';
    obj.text_s.Position = [187 170 176 19];
    obj.text_s.Text = 'Entropy [kJ/(kg K)]';

    % Create text_sR
    obj.text_sR = uieditfield(obj.Panel_parameters, 'numeric');
    obj.text_sR.ValueDisplayFormat = '%.4g';
    obj.text_sR.Editable = 'off';
    obj.text_sR.HorizontalAlignment = 'center';
    obj.text_sR.Position = [86 170 91 19];

    % Create text_sP
    obj.text_sP = uieditfield(obj.Panel_parameters, 'numeric');
    obj.text_sP.ValueDisplayFormat = '%.4g';
    obj.text_sP.Editable = 'off';
    obj.text_sP.HorizontalAlignment = 'center';
    obj.text_sP.Position = [373 170 91 19];

    % Create text_gamma
    obj.text_gamma = uilabel(obj.Panel_parameters);
    obj.text_gamma.HorizontalAlignment = 'center';
    obj.text_gamma.Position = [187 119 176 19];
    obj.text_gamma.Text = 'Adiabatic index [-]';

    % Create text_gammaR
    obj.text_gammaR = uieditfield(obj.Panel_parameters, 'numeric');
    obj.text_gammaR.ValueDisplayFormat = '%.4g';
    obj.text_gammaR.Editable = 'off';
    obj.text_gammaR.HorizontalAlignment = 'center';
    obj.text_gammaR.Position = [86 119 91 19];

    % Create text_gammaP
    obj.text_gammaP = uieditfield(obj.Panel_parameters, 'numeric');
    obj.text_gammaP.ValueDisplayFormat = '%.4g';
    obj.text_gammaP.Editable = 'off';
    obj.text_gammaP.HorizontalAlignment = 'center';
    obj.text_gammaP.Position = [373 119 91 19];

    % Create text_M
    obj.text_M = uilabel(obj.Panel_parameters);
    obj.text_M.HorizontalAlignment = 'center';
    obj.text_M.Visible = 'off';
    obj.text_M.Position = [187 16 176 19];
    obj.text_M.Text = 'Mach number [-]';

    % Create text_MP
    obj.text_MP = uieditfield(obj.Panel_parameters, 'numeric');
    obj.text_MP.ValueDisplayFormat = '%.4g';
    obj.text_MP.Editable = 'off';
    obj.text_MP.HorizontalAlignment = 'center';
    obj.text_MP.Visible = 'off';
    obj.text_MP.Position = [373 16 91 19];

    % Create text_u
    obj.text_u = uilabel(obj.Panel_parameters);
    obj.text_u.HorizontalAlignment = 'center';
    obj.text_u.Visible = 'off';
    obj.text_u.Position = [187 41 176 19];
    obj.text_u.Text = 'Flow velocity [m/s]';

    % Create text_uP
    obj.text_uP = uieditfield(obj.Panel_parameters, 'numeric');
    obj.text_uP.ValueDisplayFormat = '%.4g';
    obj.text_uP.Editable = 'off';
    obj.text_uP.HorizontalAlignment = 'center';
    obj.text_uP.Visible = 'off';
    obj.text_uP.Position = [373 41 91 19];

    % Create text_e
    obj.text_e = uilabel(obj.Panel_parameters);
    obj.text_e.HorizontalAlignment = 'center';
    obj.text_e.Position = [187 196 176 19];
    obj.text_e.Text = 'Internal energy [kJ/kg]';

    % Create text_eR
    obj.text_eR = uieditfield(obj.Panel_parameters, 'numeric');
    obj.text_eR.ValueDisplayFormat = '%.4g';
    obj.text_eR.Editable = 'off';
    obj.text_eR.HorizontalAlignment = 'center';
    obj.text_eR.Position = [86 196 91 19];

    % Create text_eP
    obj.text_eP = uieditfield(obj.Panel_parameters, 'numeric');
    obj.text_eP.ValueDisplayFormat = '%.4g';
    obj.text_eP.Editable = 'off';
    obj.text_eP.HorizontalAlignment = 'center';
    obj.text_eP.Position = [373 196 91 19];

    % Create text_sound
    obj.text_sound = uilabel(obj.Panel_parameters);
    obj.text_sound.HorizontalAlignment = 'center';
    obj.text_sound.Position = [187 67 176 19];
    obj.text_sound.Text = 'Sound speed [m/s]';

    % Create text_soundR
    obj.text_soundR = uieditfield(obj.Panel_parameters, 'numeric');
    obj.text_soundR.ValueDisplayFormat = '%.4g';
    obj.text_soundR.Editable = 'off';
    obj.text_soundR.HorizontalAlignment = 'center';
    obj.text_soundR.Position = [86 67 91 19];

    % Create text_soundP
    obj.text_soundP = uieditfield(obj.Panel_parameters, 'numeric');
    obj.text_soundP.ValueDisplayFormat = '%.4g';
    obj.text_soundP.Editable = 'off';
    obj.text_soundP.HorizontalAlignment = 'center';
    obj.text_soundP.Position = [373 67 91 19];

    % Create text_W
    obj.text_W = uilabel(obj.Panel_parameters);
    obj.text_W.HorizontalAlignment = 'center';
    obj.text_W.Position = [187 93 176 19];
    obj.text_W.Text = 'Mean Molecular Weight [g/mol]';

    % Create text_WR
    obj.text_WR = uieditfield(obj.Panel_parameters, 'numeric');
    obj.text_WR.ValueDisplayFormat = '%.4g';
    obj.text_WR.Editable = 'off';
    obj.text_WR.HorizontalAlignment = 'center';
    obj.text_WR.Position = [86 93 91 19];

    % Create text_WP
    obj.text_WP = uieditfield(obj.Panel_parameters, 'numeric');
    obj.text_WP.ValueDisplayFormat = '%.4g';
    obj.text_WP.Editable = 'off';
    obj.text_WP.HorizontalAlignment = 'center';
    obj.text_WP.Position = [373 93 91 19];

    % Create text_cp
    obj.text_cp = uilabel(obj.Panel_parameters);
    obj.text_cp.HorizontalAlignment = 'center';
    obj.text_cp.Position = [187 144 176 19];
    obj.text_cp.Text = 'Specific heat cp [kJ/(kg K)]';

    % Create text_cpR
    obj.text_cpR = uieditfield(obj.Panel_parameters, 'numeric');
    obj.text_cpR.ValueDisplayFormat = '%.4g';
    obj.text_cpR.Editable = 'off';
    obj.text_cpR.HorizontalAlignment = 'center';
    obj.text_cpR.Position = [86 144 91 19];

    % Create text_cpP
    obj.text_cpP = uieditfield(obj.Panel_parameters, 'numeric');
    obj.text_cpP.ValueDisplayFormat = '%.4g';
    obj.text_cpP.Editable = 'off';
    obj.text_cpP.HorizontalAlignment = 'center';
    obj.text_cpP.Position = [373 144 91 19];

    % Create text_h
    obj.text_h = uilabel(obj.Panel_parameters);
    obj.text_h.HorizontalAlignment = 'center';
    obj.text_h.Position = [187 221 176 19];
    obj.text_h.Text = 'Enthalpy [kJ/kg]';

    % Create text_hR
    obj.text_hR = uieditfield(obj.Panel_parameters, 'numeric');
    obj.text_hR.ValueDisplayFormat = '%.4g';
    obj.text_hR.Editable = 'off';
    obj.text_hR.HorizontalAlignment = 'center';
    obj.text_hR.Position = [86 221 91 19];

    % Create text_hP
    obj.text_hP = uieditfield(obj.Panel_parameters, 'numeric');
    obj.text_hP.ValueDisplayFormat = '%.4g';
    obj.text_hP.Editable = 'off';
    obj.text_hP.HorizontalAlignment = 'center';
    obj.text_hP.Position = [373 221 91 19];

    % Create text_r
    obj.text_r = uilabel(obj.Panel_parameters);
    obj.text_r.HorizontalAlignment = 'center';
    obj.text_r.Position = [187 247 176 19];
    obj.text_r.Text = 'Density [kg/m3]';

    % Create text_rR
    obj.text_rR = uieditfield(obj.Panel_parameters, 'numeric');
    obj.text_rR.ValueDisplayFormat = '%.4g';
    obj.text_rR.Editable = 'off';
    obj.text_rR.HorizontalAlignment = 'center';
    obj.text_rR.Position = [86 247 91 19];

    % Create text_rP
    obj.text_rP = uieditfield(obj.Panel_parameters, 'numeric');
    obj.text_rP.ValueDisplayFormat = '%.4g';
    obj.text_rP.Editable = 'off';
    obj.text_rP.HorizontalAlignment = 'center';
    obj.text_rP.Position = [373 247 91 19];

    % Create text_p
    obj.text_p = uilabel(obj.Panel_parameters);
    obj.text_p.HorizontalAlignment = 'center';
    obj.text_p.Position = [187 273 176 19];
    obj.text_p.Text = 'Pressure [bar]';

    % Create text_pP
    obj.text_pP = uieditfield(obj.Panel_parameters, 'numeric');
    obj.text_pP.Limits = [0 Inf];
    obj.text_pP.ValueDisplayFormat = '%.4g';
    obj.text_pP.Editable = 'off';
    obj.text_pP.HorizontalAlignment = 'center';
    obj.text_pP.Position = [373 273 91 19];

    % Create text_Products
    obj.text_Products = uilabel(obj.Panel_parameters);
    obj.text_Products.HorizontalAlignment = 'center';
    obj.text_Products.FontWeight = 'bold';
    obj.text_Products.Position = [373 324 91 19];
    obj.text_Products.Text = 'Products';

    % Create text_Reactans
    obj.text_Reactans = uilabel(obj.Panel_parameters);
    obj.text_Reactans.HorizontalAlignment = 'center';
    obj.text_Reactans.FontWeight = 'bold';
    obj.text_Reactans.Position = [86 321 91 22];
    obj.text_Reactans.Text = 'Reactants';

    % Create text_T
    obj.text_T = uilabel(obj.Panel_parameters);
    obj.text_T.HorizontalAlignment = 'center';
    obj.text_T.Position = [187 299 176 19];
    obj.text_T.Text = 'Temperature [K]';

    % Create text_TR
    obj.text_TR = uieditfield(obj.Panel_parameters, 'numeric');
    obj.text_TR.ValueDisplayFormat = '%.4g';
    obj.text_TR.Editable = 'off';
    obj.text_TR.HorizontalAlignment = 'center';
    obj.text_TR.Position = [86 299 91 19];

    % Create text_pR
    obj.text_pR = uieditfield(obj.Panel_parameters, 'numeric');
    obj.text_pR.ValueDisplayFormat = '%.4g';
    obj.text_pR.Editable = 'off';
    obj.text_pR.HorizontalAlignment = 'center';
    obj.text_pR.Position = [86 273 91 19];

    % Create text_TP
    obj.text_TP = uieditfield(obj.Panel_parameters, 'numeric');
    obj.text_TP.ValueDisplayFormat = '%.4g';
    obj.text_TP.Editable = 'off';
    obj.text_TP.HorizontalAlignment = 'center';
    obj.text_TP.Position = [373 299 91 19];

    % Create text_uR
    obj.text_uR = uieditfield(obj.Panel_parameters, 'numeric');
    obj.text_uR.ValueDisplayFormat = '%.4g';
    obj.text_uR.Editable = 'off';
    obj.text_uR.HorizontalAlignment = 'center';
    obj.text_uR.Visible = 'off';
    obj.text_uR.Position = [86 41 91 19];

    % Create text_MR
    obj.text_MR = uieditfield(obj.Panel_parameters, 'numeric');
    obj.text_MR.ValueDisplayFormat = '%.4g';
    obj.text_MR.Editable = 'off';
    obj.text_MR.HorizontalAlignment = 'center';
    obj.text_MR.Visible = 'off';
    obj.text_MR.Position = [86 16 91 19];

    % Create edit_phi3
    obj.edit_phi3 = uieditfield(obj.ParametersTab, 'text');
    obj.edit_phi3.Editable = 'off';
    obj.edit_phi3.HorizontalAlignment = 'center';
    obj.edit_phi3.Position = [255 409 64 19];
    obj.edit_phi3.Value = '-';

    % Create text_phi_3
    obj.text_phi_3 = uilabel(obj.ParametersTab);
    obj.text_phi_3.HorizontalAlignment = 'center';
    obj.text_phi_3.FontWeight = 'bold';
    obj.text_phi_3.Position = [270 428 35 18];
    obj.text_phi_3.Text = 'phi';
end