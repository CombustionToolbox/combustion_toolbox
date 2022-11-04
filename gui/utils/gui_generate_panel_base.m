function gui_generate_panel_base(app)
    % Routine that generates all the components to read the thermodynamic
    % properties at different stages of the rocket
    
    % Create hPhRkJkgEditFieldLabel
    app.hPhRkJkgEditFieldLabel = uilabel(app.ParametersTab);
    app.hPhRkJkgEditFieldLabel.HorizontalAlignment = 'right';
    app.hPhRkJkgEditFieldLabel.Position = [103 46 84 22];
    app.hPhRkJkgEditFieldLabel.Text = 'hP - hR [kJ/kg]';

    % Create text_q
    app.text_q = uieditfield(app.ParametersTab, 'numeric');
    app.text_q.Editable = 'off';
    app.text_q.Position = [202 46 91 19];

    % Create VolumeProductsReactantsLabel
    app.VolumeProductsReactantsLabel = uilabel(app.ParametersTab);
    app.VolumeProductsReactantsLabel.HorizontalAlignment = 'right';
    app.VolumeProductsReactantsLabel.Position = [33 21 154 22];
    app.VolumeProductsReactantsLabel.Text = 'Volume Products/Reactants';

    % Create text_vP_vR
    app.text_vP_vR = uieditfield(app.ParametersTab, 'numeric');
    app.text_vP_vR.Editable = 'off';
    app.text_vP_vR.Position = [202 21 91 19];

    % Create EpsilonmolesLabel
    app.EpsilonmolesLabel = uilabel(app.ParametersTab);
    app.EpsilonmolesLabel.HorizontalAlignment = 'right';
    app.EpsilonmolesLabel.Position = [317 46 88 22];
    app.EpsilonmolesLabel.Text = 'Epsilon (moles)';

    % Create text_error_moles
    app.text_error_moles = uieditfield(app.ParametersTab, 'numeric');
    app.text_error_moles.ValueDisplayFormat = '%11.4e';
    app.text_error_moles.ValueChangedFcn = createCallbackFcn(app, @text_error_molesValueChanged, true);
    app.text_error_moles.Editable = 'off';
    app.text_error_moles.Position = [420 46 91 19];

    % Create EpsilonmolesLabel_2
    app.EpsilonmolesLabel_2 = uilabel(app.ParametersTab);
    app.EpsilonmolesLabel_2.HorizontalAlignment = 'right';
    app.EpsilonmolesLabel_2.Position = [309 21 96 22];
    app.EpsilonmolesLabel_2.Text = 'Epsilon (method)';

    % Create text_error_problem
    app.text_error_problem = uieditfield(app.ParametersTab, 'numeric');
    app.text_error_problem.ValueDisplayFormat = '%11.4e';
    app.text_error_problem.ValueChangedFcn = createCallbackFcn(app, @text_error_problemValueChanged, true);
    app.text_error_problem.Editable = 'off';
    app.text_error_problem.Position = [420 21 91 19];

    % Create Panel_parameters
    app.Panel_parameters = uipanel(app.ParametersTab);
    app.Panel_parameters.AutoResizeChildren = 'off';
    app.Panel_parameters.BorderType = 'none';
    app.Panel_parameters.BackgroundColor = [0.9098 0.9098 0.8902];
    app.Panel_parameters.Position = [11 86 549 357];

    % Create text_s
    app.text_s = uilabel(app.Panel_parameters);
    app.text_s.HorizontalAlignment = 'center';
    app.text_s.Position = [187 170 176 19];
    app.text_s.Text = 'Entropy [kJ/(kg K)]';

    % Create text_sR
    app.text_sR = uieditfield(app.Panel_parameters, 'numeric');
    app.text_sR.ValueDisplayFormat = '%.4g';
    app.text_sR.Editable = 'off';
    app.text_sR.HorizontalAlignment = 'center';
    app.text_sR.Position = [86 170 91 19];

    % Create text_sP
    app.text_sP = uieditfield(app.Panel_parameters, 'numeric');
    app.text_sP.ValueDisplayFormat = '%.4g';
    app.text_sP.Editable = 'off';
    app.text_sP.HorizontalAlignment = 'center';
    app.text_sP.Position = [373 170 91 19];

    % Create text_gamma
    app.text_gamma = uilabel(app.Panel_parameters);
    app.text_gamma.HorizontalAlignment = 'center';
    app.text_gamma.Position = [187 119 176 19];
    app.text_gamma.Text = 'Adiabatic index [-]';

    % Create text_gammaR
    app.text_gammaR = uieditfield(app.Panel_parameters, 'numeric');
    app.text_gammaR.ValueDisplayFormat = '%.4g';
    app.text_gammaR.Editable = 'off';
    app.text_gammaR.HorizontalAlignment = 'center';
    app.text_gammaR.Position = [86 119 91 19];

    % Create text_gammaP
    app.text_gammaP = uieditfield(app.Panel_parameters, 'numeric');
    app.text_gammaP.ValueDisplayFormat = '%.4g';
    app.text_gammaP.Editable = 'off';
    app.text_gammaP.HorizontalAlignment = 'center';
    app.text_gammaP.Position = [373 119 91 19];

    % Create text_M
    app.text_M = uilabel(app.Panel_parameters);
    app.text_M.HorizontalAlignment = 'center';
    app.text_M.Visible = 'off';
    app.text_M.Position = [187 16 176 19];
    app.text_M.Text = 'Mach number [-]';

    % Create text_MP
    app.text_MP = uieditfield(app.Panel_parameters, 'numeric');
    app.text_MP.ValueDisplayFormat = '%.4g';
    app.text_MP.Editable = 'off';
    app.text_MP.HorizontalAlignment = 'center';
    app.text_MP.Visible = 'off';
    app.text_MP.Position = [373 16 91 19];

    % Create text_u
    app.text_u = uilabel(app.Panel_parameters);
    app.text_u.HorizontalAlignment = 'center';
    app.text_u.Visible = 'off';
    app.text_u.Position = [187 41 176 19];
    app.text_u.Text = 'Flow velocity [m/s]';

    % Create text_uP
    app.text_uP = uieditfield(app.Panel_parameters, 'numeric');
    app.text_uP.ValueDisplayFormat = '%.4g';
    app.text_uP.Editable = 'off';
    app.text_uP.HorizontalAlignment = 'center';
    app.text_uP.Visible = 'off';
    app.text_uP.Position = [373 41 91 19];

    % Create text_e
    app.text_e = uilabel(app.Panel_parameters);
    app.text_e.HorizontalAlignment = 'center';
    app.text_e.Position = [187 196 176 19];
    app.text_e.Text = 'Internal energy [kJ/kg]';

    % Create text_eR
    app.text_eR = uieditfield(app.Panel_parameters, 'numeric');
    app.text_eR.ValueDisplayFormat = '%.4g';
    app.text_eR.Editable = 'off';
    app.text_eR.HorizontalAlignment = 'center';
    app.text_eR.Position = [86 196 91 19];

    % Create text_eP
    app.text_eP = uieditfield(app.Panel_parameters, 'numeric');
    app.text_eP.ValueDisplayFormat = '%.4g';
    app.text_eP.Editable = 'off';
    app.text_eP.HorizontalAlignment = 'center';
    app.text_eP.Position = [373 196 91 19];

    % Create text_sound
    app.text_sound = uilabel(app.Panel_parameters);
    app.text_sound.HorizontalAlignment = 'center';
    app.text_sound.Position = [187 67 176 19];
    app.text_sound.Text = 'Sound speed [m/s]';

    % Create text_soundR
    app.text_soundR = uieditfield(app.Panel_parameters, 'numeric');
    app.text_soundR.ValueDisplayFormat = '%.4g';
    app.text_soundR.Editable = 'off';
    app.text_soundR.HorizontalAlignment = 'center';
    app.text_soundR.Position = [86 67 91 19];

    % Create text_soundP
    app.text_soundP = uieditfield(app.Panel_parameters, 'numeric');
    app.text_soundP.ValueDisplayFormat = '%.4g';
    app.text_soundP.Editable = 'off';
    app.text_soundP.HorizontalAlignment = 'center';
    app.text_soundP.Position = [373 67 91 19];

    % Create text_W
    app.text_W = uilabel(app.Panel_parameters);
    app.text_W.HorizontalAlignment = 'center';
    app.text_W.Position = [187 93 176 19];
    app.text_W.Text = 'Mean Molecular Weight [g/mol]';

    % Create text_WR
    app.text_WR = uieditfield(app.Panel_parameters, 'numeric');
    app.text_WR.ValueDisplayFormat = '%.4g';
    app.text_WR.Editable = 'off';
    app.text_WR.HorizontalAlignment = 'center';
    app.text_WR.Position = [86 93 91 19];

    % Create text_WP
    app.text_WP = uieditfield(app.Panel_parameters, 'numeric');
    app.text_WP.ValueDisplayFormat = '%.4g';
    app.text_WP.Editable = 'off';
    app.text_WP.HorizontalAlignment = 'center';
    app.text_WP.Position = [373 93 91 19];

    % Create text_cp
    app.text_cp = uilabel(app.Panel_parameters);
    app.text_cp.HorizontalAlignment = 'center';
    app.text_cp.Position = [187 144 176 19];
    app.text_cp.Text = 'Specific heat cp [kJ/(kg K)]';

    % Create text_cpR
    app.text_cpR = uieditfield(app.Panel_parameters, 'numeric');
    app.text_cpR.ValueDisplayFormat = '%.4g';
    app.text_cpR.Editable = 'off';
    app.text_cpR.HorizontalAlignment = 'center';
    app.text_cpR.Position = [86 144 91 19];

    % Create text_cpP
    app.text_cpP = uieditfield(app.Panel_parameters, 'numeric');
    app.text_cpP.ValueDisplayFormat = '%.4g';
    app.text_cpP.Editable = 'off';
    app.text_cpP.HorizontalAlignment = 'center';
    app.text_cpP.Position = [373 144 91 19];

    % Create text_h
    app.text_h = uilabel(app.Panel_parameters);
    app.text_h.HorizontalAlignment = 'center';
    app.text_h.Position = [187 221 176 19];
    app.text_h.Text = 'Enthalpy [kJ/kg]';

    % Create text_hR
    app.text_hR = uieditfield(app.Panel_parameters, 'numeric');
    app.text_hR.ValueDisplayFormat = '%.4g';
    app.text_hR.Editable = 'off';
    app.text_hR.HorizontalAlignment = 'center';
    app.text_hR.Position = [86 221 91 19];

    % Create text_hP
    app.text_hP = uieditfield(app.Panel_parameters, 'numeric');
    app.text_hP.ValueDisplayFormat = '%.4g';
    app.text_hP.Editable = 'off';
    app.text_hP.HorizontalAlignment = 'center';
    app.text_hP.Position = [373 221 91 19];

    % Create text_r
    app.text_r = uilabel(app.Panel_parameters);
    app.text_r.HorizontalAlignment = 'center';
    app.text_r.Position = [187 247 176 19];
    app.text_r.Text = 'Density [kg/m3]';

    % Create text_rR
    app.text_rR = uieditfield(app.Panel_parameters, 'numeric');
    app.text_rR.ValueDisplayFormat = '%.4g';
    app.text_rR.Editable = 'off';
    app.text_rR.HorizontalAlignment = 'center';
    app.text_rR.Position = [86 247 91 19];

    % Create text_rP
    app.text_rP = uieditfield(app.Panel_parameters, 'numeric');
    app.text_rP.ValueDisplayFormat = '%.4g';
    app.text_rP.Editable = 'off';
    app.text_rP.HorizontalAlignment = 'center';
    app.text_rP.Position = [373 247 91 19];

    % Create text_p
    app.text_p = uilabel(app.Panel_parameters);
    app.text_p.HorizontalAlignment = 'center';
    app.text_p.Position = [187 273 176 19];
    app.text_p.Text = 'Pressure [bar]';

    % Create text_pP
    app.text_pP = uieditfield(app.Panel_parameters, 'numeric');
    app.text_pP.Limits = [0 Inf];
    app.text_pP.ValueDisplayFormat = '%.4g';
    app.text_pP.Editable = 'off';
    app.text_pP.HorizontalAlignment = 'center';
    app.text_pP.Position = [373 273 91 19];

    % Create text_Products
    app.text_Products = uilabel(app.Panel_parameters);
    app.text_Products.HorizontalAlignment = 'center';
    app.text_Products.FontWeight = 'bold';
    app.text_Products.Position = [373 324 91 19];
    app.text_Products.Text = 'Products';

    % Create text_Reactans
    app.text_Reactans = uilabel(app.Panel_parameters);
    app.text_Reactans.HorizontalAlignment = 'center';
    app.text_Reactans.FontWeight = 'bold';
    app.text_Reactans.Position = [86 321 91 22];
    app.text_Reactans.Text = 'Reactants';

    % Create text_T
    app.text_T = uilabel(app.Panel_parameters);
    app.text_T.HorizontalAlignment = 'center';
    app.text_T.Position = [187 299 176 19];
    app.text_T.Text = 'Temperature [K]';

    % Create text_TR
    app.text_TR = uieditfield(app.Panel_parameters, 'numeric');
    app.text_TR.ValueDisplayFormat = '%.4g';
    app.text_TR.Editable = 'off';
    app.text_TR.HorizontalAlignment = 'center';
    app.text_TR.Position = [86 299 91 19];

    % Create text_pR
    app.text_pR = uieditfield(app.Panel_parameters, 'numeric');
    app.text_pR.ValueDisplayFormat = '%.4g';
    app.text_pR.Editable = 'off';
    app.text_pR.HorizontalAlignment = 'center';
    app.text_pR.Position = [86 273 91 19];

    % Create text_TP
    app.text_TP = uieditfield(app.Panel_parameters, 'numeric');
    app.text_TP.ValueDisplayFormat = '%.4g';
    app.text_TP.Editable = 'off';
    app.text_TP.HorizontalAlignment = 'center';
    app.text_TP.Position = [373 299 91 19];

    % Create text_uR
    app.text_uR = uieditfield(app.Panel_parameters, 'numeric');
    app.text_uR.ValueDisplayFormat = '%.4g';
    app.text_uR.Editable = 'off';
    app.text_uR.HorizontalAlignment = 'center';
    app.text_uR.Visible = 'off';
    app.text_uR.Position = [86 41 91 19];

    % Create text_MR
    app.text_MR = uieditfield(app.Panel_parameters, 'numeric');
    app.text_MR.ValueDisplayFormat = '%.4g';
    app.text_MR.Editable = 'off';
    app.text_MR.HorizontalAlignment = 'center';
    app.text_MR.Visible = 'off';
    app.text_MR.Position = [86 16 91 19];

    % Create edit_phi3
    app.edit_phi3 = uieditfield(app.ParametersTab, 'text');
    app.edit_phi3.Editable = 'off';
    app.edit_phi3.HorizontalAlignment = 'center';
    app.edit_phi3.Position = [255 409 64 19];
    app.edit_phi3.Value = '-';

    % Create text_phi_3
    app.text_phi_3 = uilabel(app.ParametersTab);
    app.text_phi_3.HorizontalAlignment = 'center';
    app.text_phi_3.FontWeight = 'bold';
    app.text_phi_3.Position = [270 428 35 18];
    app.text_phi_3.Text = 'phi';
end