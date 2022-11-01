function gui_generate_panel_rocket(obj)
    % Routine that generates all the components to read the thermodynamic
    % properties at different stages of the rocket

    % Create Panel_extra_4
    obj.Panel_extra_4 = uipanel(obj.ParametersRocketsTab);
    obj.Panel_extra_4.AutoResizeChildren = 'off';
    obj.Panel_extra_4.BorderType = 'none';
    obj.Panel_extra_4.BackgroundColor = [0.9098 0.9098 0.8902];
    obj.Panel_extra_4.Position = [457 8 100 433];

    % Create text_sP_4
    obj.text_sP_4 = uieditfield(obj.Panel_extra_4, 'numeric');
    obj.text_sP_4.ValueDisplayFormat = '%.4g';
    obj.text_sP_4.Editable = 'off';
    obj.text_sP_4.HorizontalAlignment = 'center';
    obj.text_sP_4.Position = [5 258 91 19];

    % Create text_gammaP_4
    obj.text_gammaP_4 = uieditfield(obj.Panel_extra_4, 'numeric');
    obj.text_gammaP_4.ValueDisplayFormat = '%.4g';
    obj.text_gammaP_4.Editable = 'off';
    obj.text_gammaP_4.HorizontalAlignment = 'center';
    obj.text_gammaP_4.Position = [5 207 91 19];

    % Create text_MP_4
    obj.text_MP_4 = uieditfield(obj.Panel_extra_4, 'numeric');
    obj.text_MP_4.ValueDisplayFormat = '%.4g';
    obj.text_MP_4.Editable = 'off';
    obj.text_MP_4.HorizontalAlignment = 'center';
    obj.text_MP_4.Position = [5 104 91 19];

    % Create text_uP_4
    obj.text_uP_4 = uieditfield(obj.Panel_extra_4, 'numeric');
    obj.text_uP_4.ValueDisplayFormat = '%.4g';
    obj.text_uP_4.Editable = 'off';
    obj.text_uP_4.HorizontalAlignment = 'center';
    obj.text_uP_4.Position = [5 129 91 19];

    % Create text_eP_4
    obj.text_eP_4 = uieditfield(obj.Panel_extra_4, 'numeric');
    obj.text_eP_4.ValueDisplayFormat = '%.4g';
    obj.text_eP_4.Editable = 'off';
    obj.text_eP_4.HorizontalAlignment = 'center';
    obj.text_eP_4.Position = [5 284 91 19];

    % Create text_soundP_4
    obj.text_soundP_4 = uieditfield(obj.Panel_extra_4, 'numeric');
    obj.text_soundP_4.ValueDisplayFormat = '%.4g';
    obj.text_soundP_4.Editable = 'off';
    obj.text_soundP_4.HorizontalAlignment = 'center';
    obj.text_soundP_4.Position = [5 155 91 19];

    % Create text_WP_4
    obj.text_WP_4 = uieditfield(obj.Panel_extra_4, 'numeric');
    obj.text_WP_4.ValueDisplayFormat = '%.4g';
    obj.text_WP_4.Editable = 'off';
    obj.text_WP_4.HorizontalAlignment = 'center';
    obj.text_WP_4.Position = [5 181 91 19];

    % Create text_cpP_4
    obj.text_cpP_4 = uieditfield(obj.Panel_extra_4, 'numeric');
    obj.text_cpP_4.ValueDisplayFormat = '%.4g';
    obj.text_cpP_4.Editable = 'off';
    obj.text_cpP_4.HorizontalAlignment = 'center';
    obj.text_cpP_4.Position = [5 232 91 19];

    % Create text_hP_4
    obj.text_hP_4 = uieditfield(obj.Panel_extra_4, 'numeric');
    obj.text_hP_4.ValueDisplayFormat = '%.4g';
    obj.text_hP_4.Editable = 'off';
    obj.text_hP_4.HorizontalAlignment = 'center';
    obj.text_hP_4.Position = [5 309 91 19];

    % Create text_rP_4
    obj.text_rP_4 = uieditfield(obj.Panel_extra_4, 'numeric');
    obj.text_rP_4.ValueDisplayFormat = '%.4g';
    obj.text_rP_4.Editable = 'off';
    obj.text_rP_4.HorizontalAlignment = 'center';
    obj.text_rP_4.Position = [5 335 91 19];

    % Create text_pP_4
    obj.text_pP_4 = uieditfield(obj.Panel_extra_4, 'numeric');
    obj.text_pP_4.Limits = [0 Inf];
    obj.text_pP_4.ValueDisplayFormat = '%.4g';
    obj.text_pP_4.Editable = 'off';
    obj.text_pP_4.HorizontalAlignment = 'center';
    obj.text_pP_4.Position = [5 361 91 19];

    % Create text_Products_4
    obj.text_Products_4 = uilabel(obj.Panel_extra_4);
    obj.text_Products_4.HorizontalAlignment = 'center';
    obj.text_Products_4.VerticalAlignment = 'top';
    obj.text_Products_4.FontWeight = 'bold';
    obj.text_Products_4.Position = [5 409 91 22];
    obj.text_Products_4.Text = 'Exit';

    % Create text_TP_4
    obj.text_TP_4 = uieditfield(obj.Panel_extra_4, 'numeric');
    obj.text_TP_4.ValueDisplayFormat = '%.4g';
    obj.text_TP_4.Editable = 'off';
    obj.text_TP_4.HorizontalAlignment = 'center';
    obj.text_TP_4.Position = [5 387 91 19];

    % Create text_Aratio_4
    obj.text_Aratio_4 = uieditfield(obj.Panel_extra_4, 'numeric');
    obj.text_Aratio_4.ValueDisplayFormat = '%.4g';
    obj.text_Aratio_4.Editable = 'off';
    obj.text_Aratio_4.HorizontalAlignment = 'center';
    obj.text_Aratio_4.Position = [5 79 91 19];

    % Create text_Cstar_4
    obj.text_Cstar_4 = uieditfield(obj.Panel_extra_4, 'numeric');
    obj.text_Cstar_4.ValueDisplayFormat = '%.4g';
    obj.text_Cstar_4.Editable = 'off';
    obj.text_Cstar_4.HorizontalAlignment = 'center';
    obj.text_Cstar_4.Position = [5 54 91 19];

    % Create text_Ivac_4
    obj.text_Ivac_4 = uieditfield(obj.Panel_extra_4, 'numeric');
    obj.text_Ivac_4.ValueDisplayFormat = '%.4g';
    obj.text_Ivac_4.Editable = 'off';
    obj.text_Ivac_4.HorizontalAlignment = 'center';
    obj.text_Ivac_4.Position = [5 30 91 19];

    % Create text_Isp_4
    obj.text_Isp_4 = uieditfield(obj.Panel_extra_4, 'numeric');
    obj.text_Isp_4.ValueDisplayFormat = '%.4g';
    obj.text_Isp_4.Editable = 'off';
    obj.text_Isp_4.HorizontalAlignment = 'center';
    obj.text_Isp_4.Position = [5 5 91 19];

    % Create Panel_extra_3
    obj.Panel_extra_3 = uipanel(obj.ParametersRocketsTab);
    obj.Panel_extra_3.AutoResizeChildren = 'off';
    obj.Panel_extra_3.BorderType = 'none';
    obj.Panel_extra_3.BackgroundColor = [0.9098 0.9098 0.8902];
    obj.Panel_extra_3.Position = [359 9 100 432];

    % Create text_sP_3
    obj.text_sP_3 = uieditfield(obj.Panel_extra_3, 'numeric');
    obj.text_sP_3.ValueDisplayFormat = '%.4g';
    obj.text_sP_3.Editable = 'off';
    obj.text_sP_3.HorizontalAlignment = 'center';
    obj.text_sP_3.Position = [3 257 91 19];

    % Create text_gammaP_3
    obj.text_gammaP_3 = uieditfield(obj.Panel_extra_3, 'numeric');
    obj.text_gammaP_3.ValueDisplayFormat = '%.4g';
    obj.text_gammaP_3.Editable = 'off';
    obj.text_gammaP_3.HorizontalAlignment = 'center';
    obj.text_gammaP_3.Position = [3 206 91 19];

    % Create text_MP_3
    obj.text_MP_3 = uieditfield(obj.Panel_extra_3, 'numeric');
    obj.text_MP_3.ValueDisplayFormat = '%.4g';
    obj.text_MP_3.Editable = 'off';
    obj.text_MP_3.HorizontalAlignment = 'center';
    obj.text_MP_3.Position = [3 103 91 19];

    % Create text_uP_3
    obj.text_uP_3 = uieditfield(obj.Panel_extra_3, 'numeric');
    obj.text_uP_3.ValueDisplayFormat = '%.4g';
    obj.text_uP_3.Editable = 'off';
    obj.text_uP_3.HorizontalAlignment = 'center';
    obj.text_uP_3.Position = [3 128 91 19];

    % Create text_eP_3
    obj.text_eP_3 = uieditfield(obj.Panel_extra_3, 'numeric');
    obj.text_eP_3.ValueDisplayFormat = '%.4g';
    obj.text_eP_3.Editable = 'off';
    obj.text_eP_3.HorizontalAlignment = 'center';
    obj.text_eP_3.Position = [3 283 91 19];

    % Create text_soundP_3
    obj.text_soundP_3 = uieditfield(obj.Panel_extra_3, 'numeric');
    obj.text_soundP_3.ValueDisplayFormat = '%.4g';
    obj.text_soundP_3.Editable = 'off';
    obj.text_soundP_3.HorizontalAlignment = 'center';
    obj.text_soundP_3.Position = [3 154 91 19];

    % Create text_WP_3
    obj.text_WP_3 = uieditfield(obj.Panel_extra_3, 'numeric');
    obj.text_WP_3.ValueDisplayFormat = '%.4g';
    obj.text_WP_3.Editable = 'off';
    obj.text_WP_3.HorizontalAlignment = 'center';
    obj.text_WP_3.Position = [3 180 91 19];

    % Create text_cpP_3
    obj.text_cpP_3 = uieditfield(obj.Panel_extra_3, 'numeric');
    obj.text_cpP_3.ValueDisplayFormat = '%.4g';
    obj.text_cpP_3.Editable = 'off';
    obj.text_cpP_3.HorizontalAlignment = 'center';
    obj.text_cpP_3.Position = [3 231 91 19];

    % Create text_hP_3
    obj.text_hP_3 = uieditfield(obj.Panel_extra_3, 'numeric');
    obj.text_hP_3.ValueDisplayFormat = '%.4g';
    obj.text_hP_3.Editable = 'off';
    obj.text_hP_3.HorizontalAlignment = 'center';
    obj.text_hP_3.Position = [3 308 91 19];

    % Create text_rP_3
    obj.text_rP_3 = uieditfield(obj.Panel_extra_3, 'numeric');
    obj.text_rP_3.ValueDisplayFormat = '%.4g';
    obj.text_rP_3.Editable = 'off';
    obj.text_rP_3.HorizontalAlignment = 'center';
    obj.text_rP_3.Position = [3 334 91 19];

    % Create text_pP_3
    obj.text_pP_3 = uieditfield(obj.Panel_extra_3, 'numeric');
    obj.text_pP_3.Limits = [0 Inf];
    obj.text_pP_3.ValueDisplayFormat = '%.4g';
    obj.text_pP_3.Editable = 'off';
    obj.text_pP_3.HorizontalAlignment = 'center';
    obj.text_pP_3.Position = [3 360 91 19];

    % Create text_Products_3
    obj.text_Products_3 = uilabel(obj.Panel_extra_3);
    obj.text_Products_3.HorizontalAlignment = 'center';
    obj.text_Products_3.VerticalAlignment = 'top';
    obj.text_Products_3.FontWeight = 'bold';
    obj.text_Products_3.Position = [3 408 91 22];
    obj.text_Products_3.Text = 'Throat';

    % Create text_TP_3
    obj.text_TP_3 = uieditfield(obj.Panel_extra_3, 'numeric');
    obj.text_TP_3.ValueDisplayFormat = '%.4g';
    obj.text_TP_3.Editable = 'off';
    obj.text_TP_3.HorizontalAlignment = 'center';
    obj.text_TP_3.Position = [3 386 91 19];

    % Create text_Aratio_3
    obj.text_Aratio_3 = uieditfield(obj.Panel_extra_3, 'numeric');
    obj.text_Aratio_3.ValueDisplayFormat = '%.4g';
    obj.text_Aratio_3.Editable = 'off';
    obj.text_Aratio_3.HorizontalAlignment = 'center';
    obj.text_Aratio_3.Position = [3 78 91 19];

    % Create text_Cstar_3
    obj.text_Cstar_3 = uieditfield(obj.Panel_extra_3, 'numeric');
    obj.text_Cstar_3.ValueDisplayFormat = '%.4g';
    obj.text_Cstar_3.Editable = 'off';
    obj.text_Cstar_3.HorizontalAlignment = 'center';
    obj.text_Cstar_3.Position = [3 53 91 19];

    % Create text_Ivac_3
    obj.text_Ivac_3 = uieditfield(obj.Panel_extra_3, 'numeric');
    obj.text_Ivac_3.ValueDisplayFormat = '%.4g';
    obj.text_Ivac_3.Editable = 'off';
    obj.text_Ivac_3.HorizontalAlignment = 'center';
    obj.text_Ivac_3.Position = [3 29 91 19];

    % Create text_Isp_3
    obj.text_Isp_3 = uieditfield(obj.Panel_extra_3, 'numeric');
    obj.text_Isp_3.ValueDisplayFormat = '%.4g';
    obj.text_Isp_3.Editable = 'off';
    obj.text_Isp_3.HorizontalAlignment = 'center';
    obj.text_Isp_3.Position = [3 4 91 19];

    % Create Panel_extra_5
    obj.Panel_extra_5 = uipanel(obj.ParametersRocketsTab);
    obj.Panel_extra_5.AutoResizeChildren = 'off';
    obj.Panel_extra_5.BorderType = 'none';
    obj.Panel_extra_5.Visible = 'off';
    obj.Panel_extra_5.BackgroundColor = [0.9098 0.9098 0.8902];
    obj.Panel_extra_5.Position = [557 8 100 433];

    % Create text_sP_5
    obj.text_sP_5 = uieditfield(obj.Panel_extra_5, 'numeric');
    obj.text_sP_5.ValueDisplayFormat = '%.4g';
    obj.text_sP_5.Editable = 'off';
    obj.text_sP_5.HorizontalAlignment = 'center';
    obj.text_sP_5.Position = [5 258 91 19];

    % Create text_gammaP_5
    obj.text_gammaP_5 = uieditfield(obj.Panel_extra_5, 'numeric');
    obj.text_gammaP_5.ValueDisplayFormat = '%.4g';
    obj.text_gammaP_5.Editable = 'off';
    obj.text_gammaP_5.HorizontalAlignment = 'center';
    obj.text_gammaP_5.Position = [5 207 91 19];

    % Create text_MP_5
    obj.text_MP_5 = uieditfield(obj.Panel_extra_5, 'numeric');
    obj.text_MP_5.ValueDisplayFormat = '%.4g';
    obj.text_MP_5.Editable = 'off';
    obj.text_MP_5.HorizontalAlignment = 'center';
    obj.text_MP_5.Position = [5 104 91 19];

    % Create text_uP_5
    obj.text_uP_5 = uieditfield(obj.Panel_extra_5, 'numeric');
    obj.text_uP_5.ValueDisplayFormat = '%.4g';
    obj.text_uP_5.Editable = 'off';
    obj.text_uP_5.HorizontalAlignment = 'center';
    obj.text_uP_5.Position = [5 129 91 19];

    % Create text_eP_5
    obj.text_eP_5 = uieditfield(obj.Panel_extra_5, 'numeric');
    obj.text_eP_5.ValueDisplayFormat = '%.4g';
    obj.text_eP_5.Editable = 'off';
    obj.text_eP_5.HorizontalAlignment = 'center';
    obj.text_eP_5.Position = [5 284 91 19];

    % Create text_soundP_5
    obj.text_soundP_5 = uieditfield(obj.Panel_extra_5, 'numeric');
    obj.text_soundP_5.ValueDisplayFormat = '%.4g';
    obj.text_soundP_5.Editable = 'off';
    obj.text_soundP_5.HorizontalAlignment = 'center';
    obj.text_soundP_5.Position = [5 155 91 19];

    % Create text_WP_5
    obj.text_WP_5 = uieditfield(obj.Panel_extra_5, 'numeric');
    obj.text_WP_5.ValueDisplayFormat = '%.4g';
    obj.text_WP_5.Editable = 'off';
    obj.text_WP_5.HorizontalAlignment = 'center';
    obj.text_WP_5.Position = [5 181 91 19];

    % Create text_cpP_5
    obj.text_cpP_5 = uieditfield(obj.Panel_extra_5, 'numeric');
    obj.text_cpP_5.ValueDisplayFormat = '%.4g';
    obj.text_cpP_5.Editable = 'off';
    obj.text_cpP_5.HorizontalAlignment = 'center';
    obj.text_cpP_5.Position = [5 232 91 19];

    % Create text_hP_5
    obj.text_hP_5 = uieditfield(obj.Panel_extra_5, 'numeric');
    obj.text_hP_5.ValueDisplayFormat = '%.4g';
    obj.text_hP_5.Editable = 'off';
    obj.text_hP_5.HorizontalAlignment = 'center';
    obj.text_hP_5.Position = [5 309 91 19];

    % Create text_rP_5
    obj.text_rP_5 = uieditfield(obj.Panel_extra_5, 'numeric');
    obj.text_rP_5.ValueDisplayFormat = '%.4g';
    obj.text_rP_5.Editable = 'off';
    obj.text_rP_5.HorizontalAlignment = 'center';
    obj.text_rP_5.Position = [5 335 91 19];

    % Create text_pP_5
    obj.text_pP_5 = uieditfield(obj.Panel_extra_5, 'numeric');
    obj.text_pP_5.Limits = [0 Inf];
    obj.text_pP_5.ValueDisplayFormat = '%.4g';
    obj.text_pP_5.Editable = 'off';
    obj.text_pP_5.HorizontalAlignment = 'center';
    obj.text_pP_5.Position = [5 361 91 19];

    % Create text_Products_5
    obj.text_Products_5 = uilabel(obj.Panel_extra_5);
    obj.text_Products_5.HorizontalAlignment = 'center';
    obj.text_Products_5.VerticalAlignment = 'top';
    obj.text_Products_5.FontWeight = 'bold';
    obj.text_Products_5.Position = [5 409 91 22];
    obj.text_Products_5.Text = 'Exit';

    % Create text_TP_5
    obj.text_TP_5 = uieditfield(obj.Panel_extra_5, 'numeric');
    obj.text_TP_5.ValueDisplayFormat = '%.4g';
    obj.text_TP_5.Editable = 'off';
    obj.text_TP_5.HorizontalAlignment = 'center';
    obj.text_TP_5.Position = [5 387 91 19];

    % Create text_Aratio_5
    obj.text_Aratio_5 = uieditfield(obj.Panel_extra_5, 'numeric');
    obj.text_Aratio_5.ValueDisplayFormat = '%.4g';
    obj.text_Aratio_5.Editable = 'off';
    obj.text_Aratio_5.HorizontalAlignment = 'center';
    obj.text_Aratio_5.Position = [5 79 91 19];

    % Create text_Cstar_5
    obj.text_Cstar_5 = uieditfield(obj.Panel_extra_5, 'numeric');
    obj.text_Cstar_5.ValueDisplayFormat = '%.4g';
    obj.text_Cstar_5.Editable = 'off';
    obj.text_Cstar_5.HorizontalAlignment = 'center';
    obj.text_Cstar_5.Position = [5 54 91 19];

    % Create text_Ivac_5
    obj.text_Ivac_5 = uieditfield(obj.Panel_extra_5, 'numeric');
    obj.text_Ivac_5.ValueDisplayFormat = '%.4g';
    obj.text_Ivac_5.Editable = 'off';
    obj.text_Ivac_5.HorizontalAlignment = 'center';
    obj.text_Ivac_5.Position = [5 30 91 19];

    % Create text_Isp_5
    obj.text_Isp_5 = uieditfield(obj.Panel_extra_5, 'numeric');
    obj.text_Isp_5.ValueDisplayFormat = '%.4g';
    obj.text_Isp_5.Editable = 'off';
    obj.text_Isp_5.HorizontalAlignment = 'center';
    obj.text_Isp_5.Position = [5 5 91 19];

    % Create Panel_extra_0
    obj.Panel_extra_0 = uipanel(obj.ParametersRocketsTab);
    obj.Panel_extra_0.AutoResizeChildren = 'off';
    obj.Panel_extra_0.BorderType = 'none';
    obj.Panel_extra_0.BackgroundColor = [0.9098 0.9098 0.8902];
    obj.Panel_extra_0.Position = [13 9 342 436];

    % Create text_s_2
    obj.text_s_2 = uilabel(obj.Panel_extra_0);
    obj.text_s_2.HorizontalAlignment = 'center';
    obj.text_s_2.Position = [104 257 136 19];
    obj.text_s_2.Text = 'Entropy [kJ/(kg K)]';

    % Create text_sR_2
    obj.text_sR_2 = uieditfield(obj.Panel_extra_0, 'numeric');
    obj.text_sR_2.ValueDisplayFormat = '%.4g';
    obj.text_sR_2.Editable = 'off';
    obj.text_sR_2.HorizontalAlignment = 'center';
    obj.text_sR_2.Position = [6 257 91 19];

    % Create text_sP_2
    obj.text_sP_2 = uieditfield(obj.Panel_extra_0, 'numeric');
    obj.text_sP_2.ValueDisplayFormat = '%.4g';
    obj.text_sP_2.Editable = 'off';
    obj.text_sP_2.HorizontalAlignment = 'center';
    obj.text_sP_2.Position = [248 257 91 19];

    % Create text_gamma_2
    obj.text_gamma_2 = uilabel(obj.Panel_extra_0);
    obj.text_gamma_2.HorizontalAlignment = 'center';
    obj.text_gamma_2.Position = [104 206 136 19];
    obj.text_gamma_2.Text = 'Adibatic index [-]';

    % Create text_gammaR_2
    obj.text_gammaR_2 = uieditfield(obj.Panel_extra_0, 'numeric');
    obj.text_gammaR_2.ValueDisplayFormat = '%.4g';
    obj.text_gammaR_2.Editable = 'off';
    obj.text_gammaR_2.HorizontalAlignment = 'center';
    obj.text_gammaR_2.Position = [6 206 91 19];

    % Create text_gammaP_2
    obj.text_gammaP_2 = uieditfield(obj.Panel_extra_0, 'numeric');
    obj.text_gammaP_2.ValueDisplayFormat = '%.4g';
    obj.text_gammaP_2.Editable = 'off';
    obj.text_gammaP_2.HorizontalAlignment = 'center';
    obj.text_gammaP_2.Position = [248 206 91 19];

    % Create text_M_2
    obj.text_M_2 = uilabel(obj.Panel_extra_0);
    obj.text_M_2.HorizontalAlignment = 'center';
    obj.text_M_2.Position = [104 103 136 19];
    obj.text_M_2.Text = 'Mach number [-]';

    % Create text_MP_2
    obj.text_MP_2 = uieditfield(obj.Panel_extra_0, 'numeric');
    obj.text_MP_2.ValueDisplayFormat = '%.4g';
    obj.text_MP_2.Editable = 'off';
    obj.text_MP_2.HorizontalAlignment = 'center';
    obj.text_MP_2.Position = [248 103 91 19];

    % Create text_u_2
    obj.text_u_2 = uilabel(obj.Panel_extra_0);
    obj.text_u_2.HorizontalAlignment = 'center';
    obj.text_u_2.Position = [104 128 136 19];
    obj.text_u_2.Text = 'Flow velocity [m/s]';

    % Create text_uP_2
    obj.text_uP_2 = uieditfield(obj.Panel_extra_0, 'numeric');
    obj.text_uP_2.ValueDisplayFormat = '%.4g';
    obj.text_uP_2.Editable = 'off';
    obj.text_uP_2.HorizontalAlignment = 'center';
    obj.text_uP_2.Position = [248 128 91 19];

    % Create text_e_2
    obj.text_e_2 = uilabel(obj.Panel_extra_0);
    obj.text_e_2.HorizontalAlignment = 'center';
    obj.text_e_2.Position = [104 283 136 19];
    obj.text_e_2.Text = 'Internal energy [kJ/kg]';

    % Create text_eR_2
    obj.text_eR_2 = uieditfield(obj.Panel_extra_0, 'numeric');
    obj.text_eR_2.ValueDisplayFormat = '%.4g';
    obj.text_eR_2.Editable = 'off';
    obj.text_eR_2.HorizontalAlignment = 'center';
    obj.text_eR_2.Position = [6 283 91 19];

    % Create text_eP_2
    obj.text_eP_2 = uieditfield(obj.Panel_extra_0, 'numeric');
    obj.text_eP_2.ValueDisplayFormat = '%.4g';
    obj.text_eP_2.Editable = 'off';
    obj.text_eP_2.HorizontalAlignment = 'center';
    obj.text_eP_2.Position = [248 283 91 19];

    % Create text_sound_2
    obj.text_sound_2 = uilabel(obj.Panel_extra_0);
    obj.text_sound_2.HorizontalAlignment = 'center';
    obj.text_sound_2.Position = [104 154 136 19];
    obj.text_sound_2.Text = 'Sound speed [m/s]';

    % Create text_soundR_2
    obj.text_soundR_2 = uieditfield(obj.Panel_extra_0, 'numeric');
    obj.text_soundR_2.ValueDisplayFormat = '%.4g';
    obj.text_soundR_2.Editable = 'off';
    obj.text_soundR_2.HorizontalAlignment = 'center';
    obj.text_soundR_2.Position = [6 154 91 19];

    % Create text_soundP_2
    obj.text_soundP_2 = uieditfield(obj.Panel_extra_0, 'numeric');
    obj.text_soundP_2.ValueDisplayFormat = '%.4g';
    obj.text_soundP_2.Editable = 'off';
    obj.text_soundP_2.HorizontalAlignment = 'center';
    obj.text_soundP_2.Position = [248 154 91 19];

    % Create text_W_2
    obj.text_W_2 = uilabel(obj.Panel_extra_0);
    obj.text_W_2.HorizontalAlignment = 'center';
    obj.text_W_2.Position = [104 180 136 19];
    obj.text_W_2.Text = 'Molecular Weight [g/mol]';

    % Create text_WR_2
    obj.text_WR_2 = uieditfield(obj.Panel_extra_0, 'numeric');
    obj.text_WR_2.ValueDisplayFormat = '%.4g';
    obj.text_WR_2.Editable = 'off';
    obj.text_WR_2.HorizontalAlignment = 'center';
    obj.text_WR_2.Position = [6 180 91 19];

    % Create text_WP_2
    obj.text_WP_2 = uieditfield(obj.Panel_extra_0, 'numeric');
    obj.text_WP_2.ValueDisplayFormat = '%.4g';
    obj.text_WP_2.Editable = 'off';
    obj.text_WP_2.HorizontalAlignment = 'center';
    obj.text_WP_2.Position = [248 180 91 19];

    % Create text_cp_2
    obj.text_cp_2 = uilabel(obj.Panel_extra_0);
    obj.text_cp_2.HorizontalAlignment = 'center';
    obj.text_cp_2.Position = [99 231 147 19];
    obj.text_cp_2.Text = 'Specific heat cp [kJ/(kg K)]';

    % Create text_cpR_2
    obj.text_cpR_2 = uieditfield(obj.Panel_extra_0, 'numeric');
    obj.text_cpR_2.ValueDisplayFormat = '%.4g';
    obj.text_cpR_2.Editable = 'off';
    obj.text_cpR_2.HorizontalAlignment = 'center';
    obj.text_cpR_2.Position = [6 231 91 19];

    % Create text_cpP_2
    obj.text_cpP_2 = uieditfield(obj.Panel_extra_0, 'numeric');
    obj.text_cpP_2.ValueDisplayFormat = '%.4g';
    obj.text_cpP_2.Editable = 'off';
    obj.text_cpP_2.HorizontalAlignment = 'center';
    obj.text_cpP_2.Position = [248 231 91 19];

    % Create text_h_2
    obj.text_h_2 = uilabel(obj.Panel_extra_0);
    obj.text_h_2.HorizontalAlignment = 'center';
    obj.text_h_2.Position = [104 308 136 19];
    obj.text_h_2.Text = 'Enthalpy [kJ/kg]';

    % Create text_hR_2
    obj.text_hR_2 = uieditfield(obj.Panel_extra_0, 'numeric');
    obj.text_hR_2.ValueDisplayFormat = '%.4g';
    obj.text_hR_2.Editable = 'off';
    obj.text_hR_2.HorizontalAlignment = 'center';
    obj.text_hR_2.Position = [6 308 91 19];

    % Create text_hP_2
    obj.text_hP_2 = uieditfield(obj.Panel_extra_0, 'numeric');
    obj.text_hP_2.ValueDisplayFormat = '%.4g';
    obj.text_hP_2.Editable = 'off';
    obj.text_hP_2.HorizontalAlignment = 'center';
    obj.text_hP_2.Position = [248 308 91 19];

    % Create text_r_2
    obj.text_r_2 = uilabel(obj.Panel_extra_0);
    obj.text_r_2.HorizontalAlignment = 'center';
    obj.text_r_2.Position = [104 334 136 19];
    obj.text_r_2.Text = 'Density [kg/m3]';

    % Create text_rR_2
    obj.text_rR_2 = uieditfield(obj.Panel_extra_0, 'numeric');
    obj.text_rR_2.ValueDisplayFormat = '%.4g';
    obj.text_rR_2.Editable = 'off';
    obj.text_rR_2.HorizontalAlignment = 'center';
    obj.text_rR_2.Position = [6 334 91 19];

    % Create text_rP_2
    obj.text_rP_2 = uieditfield(obj.Panel_extra_0, 'numeric');
    obj.text_rP_2.ValueDisplayFormat = '%.4g';
    obj.text_rP_2.Editable = 'off';
    obj.text_rP_2.HorizontalAlignment = 'center';
    obj.text_rP_2.Position = [248 334 91 19];

    % Create text_p_2
    obj.text_p_2 = uilabel(obj.Panel_extra_0);
    obj.text_p_2.HorizontalAlignment = 'center';
    obj.text_p_2.Position = [104 360 136 19];
    obj.text_p_2.Text = 'Pressure [bar]';

    % Create text_pP_2
    obj.text_pP_2 = uieditfield(obj.Panel_extra_0, 'numeric');
    obj.text_pP_2.Limits = [0 Inf];
    obj.text_pP_2.ValueDisplayFormat = '%.4g';
    obj.text_pP_2.Editable = 'off';
    obj.text_pP_2.HorizontalAlignment = 'center';
    obj.text_pP_2.Position = [248 360 91 19];

    % Create text_Products_2
    obj.text_Products_2 = uilabel(obj.Panel_extra_0);
    obj.text_Products_2.HorizontalAlignment = 'center';
    obj.text_Products_2.VerticalAlignment = 'top';
    obj.text_Products_2.FontWeight = 'bold';
    obj.text_Products_2.Position = [246 408 96 22];
    obj.text_Products_2.Text = 'Outlet Chamber';

    % Create text_Reactans_2
    obj.text_Reactans_2 = uilabel(obj.Panel_extra_0);
    obj.text_Reactans_2.HorizontalAlignment = 'center';
    obj.text_Reactans_2.VerticalAlignment = 'top';
    obj.text_Reactans_2.FontWeight = 'bold';
    obj.text_Reactans_2.Position = [6 411 91 19];
    obj.text_Reactans_2.Text = 'Reactants';

    % Create text_T_2
    obj.text_T_2 = uilabel(obj.Panel_extra_0);
    obj.text_T_2.HorizontalAlignment = 'center';
    obj.text_T_2.Position = [104 386 136 19];
    obj.text_T_2.Text = 'Temperature [K]';

    % Create text_TR_2
    obj.text_TR_2 = uieditfield(obj.Panel_extra_0, 'numeric');
    obj.text_TR_2.ValueDisplayFormat = '%.4g';
    obj.text_TR_2.Editable = 'off';
    obj.text_TR_2.HorizontalAlignment = 'center';
    obj.text_TR_2.Position = [6 386 91 19];

    % Create text_pR_2
    obj.text_pR_2 = uieditfield(obj.Panel_extra_0, 'numeric');
    obj.text_pR_2.ValueDisplayFormat = '%.4g';
    obj.text_pR_2.Editable = 'off';
    obj.text_pR_2.HorizontalAlignment = 'center';
    obj.text_pR_2.Position = [6 360 91 19];

    % Create text_TP_2
    obj.text_TP_2 = uieditfield(obj.Panel_extra_0, 'numeric');
    obj.text_TP_2.ValueDisplayFormat = '%.4g';
    obj.text_TP_2.Editable = 'off';
    obj.text_TP_2.HorizontalAlignment = 'center';
    obj.text_TP_2.Position = [248 386 91 19];

    % Create text_uR_2
    obj.text_uR_2 = uieditfield(obj.Panel_extra_0, 'numeric');
    obj.text_uR_2.ValueDisplayFormat = '%.4g';
    obj.text_uR_2.Editable = 'off';
    obj.text_uR_2.HorizontalAlignment = 'center';
    obj.text_uR_2.Position = [6 128 91 19];

    % Create text_MR_2
    obj.text_MR_2 = uieditfield(obj.Panel_extra_0, 'numeric');
    obj.text_MR_2.ValueDisplayFormat = '%.4g';
    obj.text_MR_2.Editable = 'off';
    obj.text_MR_2.HorizontalAlignment = 'center';
    obj.text_MR_2.Position = [6 103 91 19];

    % Create text_Aratio
    obj.text_Aratio = uilabel(obj.Panel_extra_0);
    obj.text_Aratio.HorizontalAlignment = 'center';
    obj.text_Aratio.Position = [104 78 136 19];
    obj.text_Aratio.Text = 'Area ratio A/A_t [-]';

    % Create text_Aratio_2
    obj.text_Aratio_2 = uieditfield(obj.Panel_extra_0, 'numeric');
    obj.text_Aratio_2.ValueDisplayFormat = '%.4g';
    obj.text_Aratio_2.Editable = 'off';
    obj.text_Aratio_2.HorizontalAlignment = 'center';
    obj.text_Aratio_2.Position = [248 78 91 19];

    % Create text_Aratio_1
    obj.text_Aratio_1 = uieditfield(obj.Panel_extra_0, 'numeric');
    obj.text_Aratio_1.ValueDisplayFormat = '%.4g';
    obj.text_Aratio_1.Editable = 'off';
    obj.text_Aratio_1.HorizontalAlignment = 'center';
    obj.text_Aratio_1.Position = [6 78 91 19];

    % Create text_Cstar
    obj.text_Cstar = uilabel(obj.Panel_extra_0);
    obj.text_Cstar.HorizontalAlignment = 'center';
    obj.text_Cstar.Position = [112 53 120 19];
    obj.text_Cstar.Text = 'Charac. velocity [m/s]';

    % Create text_Cstar_2
    obj.text_Cstar_2 = uieditfield(obj.Panel_extra_0, 'numeric');
    obj.text_Cstar_2.ValueDisplayFormat = '%.4g';
    obj.text_Cstar_2.Editable = 'off';
    obj.text_Cstar_2.HorizontalAlignment = 'center';
    obj.text_Cstar_2.Position = [248 53 91 19];

    % Create text_Cstar_1
    obj.text_Cstar_1 = uieditfield(obj.Panel_extra_0, 'numeric');
    obj.text_Cstar_1.ValueDisplayFormat = '%.4g';
    obj.text_Cstar_1.Editable = 'off';
    obj.text_Cstar_1.HorizontalAlignment = 'center';
    obj.text_Cstar_1.Position = [6 53 91 19];

    % Create text_Ivac
    obj.text_Ivac = uilabel(obj.Panel_extra_0);
    obj.text_Ivac.HorizontalAlignment = 'center';
    obj.text_Ivac.Position = [107 29 130 19];
    obj.text_Ivac.Text = 'Specific impulse vac [s]';

    % Create text_Ivac_2
    obj.text_Ivac_2 = uieditfield(obj.Panel_extra_0, 'numeric');
    obj.text_Ivac_2.ValueDisplayFormat = '%.4g';
    obj.text_Ivac_2.Editable = 'off';
    obj.text_Ivac_2.HorizontalAlignment = 'center';
    obj.text_Ivac_2.Position = [248 29 91 19];

    % Create text_Ivac_1
    obj.text_Ivac_1 = uieditfield(obj.Panel_extra_0, 'numeric');
    obj.text_Ivac_1.ValueDisplayFormat = '%.4g';
    obj.text_Ivac_1.Editable = 'off';
    obj.text_Ivac_1.HorizontalAlignment = 'center';
    obj.text_Ivac_1.Position = [6 29 91 19];

    % Create text_Isp
    obj.text_Isp = uilabel(obj.Panel_extra_0);
    obj.text_Isp.HorizontalAlignment = 'center';
    obj.text_Isp.Position = [118 4 108 19];
    obj.text_Isp.Text = 'Specific impulse [s]';

    % Create text_Isp_2
    obj.text_Isp_2 = uieditfield(obj.Panel_extra_0, 'numeric');
    obj.text_Isp_2.ValueDisplayFormat = '%.4g';
    obj.text_Isp_2.Editable = 'off';
    obj.text_Isp_2.HorizontalAlignment = 'center';
    obj.text_Isp_2.Position = [248 4 91 19];

    % Create text_Isp_1
    obj.text_Isp_1 = uieditfield(obj.Panel_extra_0, 'numeric');
    obj.text_Isp_1.ValueDisplayFormat = '%.4g';
    obj.text_Isp_1.Editable = 'off';
    obj.text_Isp_1.HorizontalAlignment = 'center';
    obj.text_Isp_1.Position = [6 4 91 19];
end