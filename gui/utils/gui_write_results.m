function gui_write_results(app, results, i, varargin)
    % Update GUI with the results of the ith case solved
    
    % Default values
    FLAG_REACTANTS = false;
    % Check varargin values
    for j = 1:nargin - 4
        switch lower(varargin{j})
            case 'flag_reactants'
                FLAG_REACTANTS = varargin{j + 1};
        end
    end
    % Update properties
    update_properties(app, results, i);
    % Update mixture (tables)
    update_mixtures(app, results, i, FLAG_REACTANTS);
    % Update equivalence ratio
    update_equivalence_ratio(app, results, i);
end

% SUB-PASS FUNCTIONS
function update_properties(app, results, i)
    mix1 = results(i).mix1;
    mix2 = results(i).mix2;
    
    update_properties_common(app, mix1, '_1');
    update_properties_common(app, mix2, '_2');

    if strcmpi(results(i).ProblemType, 'TP')
        app.text_error_problem.Value = mix2.error_moles;
    else
        app.text_error_problem.Value = mix2.error_problem;
    end
    if contains(results(i).ProblemType, 'SHOCK', 'IgnoreCase', true) || contains(results(i).ProblemType, 'DET', 'IgnoreCase', true) || contains(results(i).ProblemType, 'ROCKET', 'IgnoreCase', true)
        update_properties_velocities(app, mix1, '_1')
        update_properties_velocities(app, mix2, '_2')
        
        if contains(results(i).ProblemType, 'OBLIQUE', 'IgnoreCase', true)
            update_properties_oblique(app, mix2, '_2');
        end

        if contains(results(i).ProblemType, 'ROCKET', 'IgnoreCase', true)
            app.text_Aratio_2.Value = mix2.Aratio;
            
            if isfield(results(i), 'mix2_c')
                % Combustor end (c)
                mix3 = results(i).mix2_c;
                FLAG_IAC = false;
            else
                % Throat
                mix3 = results(i).mix3;
                FLAG_IAC = true;
            end
            
            update_properties_common(app, mix3, '_3');
            update_properties_velocities(app, mix3, '_3')
            update_properties_rocket(app, mix3, '_3')
            
            if ~FLAG_IAC || (~isempty(results(i).strP) && FLAG_IAC)
                
                if ~FLAG_IAC
                    % Throat
                    mix4 = results(i).mix3;
                else
                    % Exit
                    mix4 = results(i).strP;
                end

                % Update values
                update_properties_common(app, mix4, '_4');
                update_properties_velocities(app, mix4, '_4')
                update_properties_rocket(app, mix4, '_4')

                % Get max error method
                app.text_error_problem.Value = max([mix2.error_problem, mix3.error_problem, mix4.error_problem]);

            else
                % Exit
                app.text_T_4.Value = 0;
                app.text_p_4.Value = 0;
                app.text_r_4.Value = 0;
                app.text_h_4.Value = 0;
                app.text_e_4.Value = 0;
                app.text_cp_4.Value = 0;
                app.text_s_4.Value = 0;
                app.text_gamma_4.Value = 0;
                app.text_W_4.Value = 0;
                app.text_sound_4.Value = 0;
                app.text_u_4.Value = 0;
                app.text_M_4.Value = 0;
                app.text_Aratio_4.Value = 0;
                app.text_Cstar_4.Value = 0;
                app.text_Ivac_4.Value = 0;
                app.text_Isp_4.Value = 0;
            end
            
            if FLAG_IAC
                return
            end

            if ~isempty(results(i).strP)
                % Exit
                mix5 = results(i).strP;
                
                % Update values
                update_properties_common(app, mix5, '_5');
                update_properties_velocities(app, mix5, '_5')
                update_properties_rocket(app, mix5, '_5')

                % Get max error method
                app.text_error_problem.Value = max([mix2.error_problem, mix3.error_problem, mix4.error_problem, mix5.error_problem]);
            
            else
                app.text_T_5.Value = 0;
                app.text_p_5.Value = 0;
                app.text_r_5.Value = 0;
                app.text_h_5.Value = 0;
                app.text_e_5.Value = 0;
                app.text_cp_5.Value = 0;
                app.text_s_5.Value = 0;
                app.text_gamma_5.Value = 0;
                app.text_W_5.Value = 0;
                app.text_sound_5.Value = 0;
                app.text_u_5.Value = 0;
                app.text_M_5.Value = 0;
                app.text_Aratio_5.Value = 0;
                app.text_Cstar_5.Value = 0;
                app.text_Ivac_5.Value = 0;
                app.text_Isp_5.Value = 0;

                % Get max error method
                app.text_error_problem.Value = max([mix2.error_problem, mix3.error_problem, mix4.error_problem]);
            end

        end
    end
end

function update_properties_common(app, mix, suffix)
    % Update common mixture properties
    app.(['text_T', suffix]).Value = temperature(mix);
    app.(['text_p', suffix]).Value = pressure(mix);
    app.(['text_r', suffix]).Value = density(mix);
    app.(['text_h', suffix]).Value = enthalpy_mass(mix);
    app.(['text_e', suffix]).Value = intEnergy_mass(mix);
    app.(['text_cp', suffix]).Value = cp_mass(mix);
    try
        app.(['text_s', suffix]).Value = entropy_mass(mix);
    catch
        app.(['text_s', suffix]).Value = Inf;
    end
    app.(['text_gamma', suffix]).Value = adiabaticIndex(mix);
    app.(['text_W', suffix]).Value = MolecularWeight(mix);
    app.(['text_sound', suffix]).Value = soundspeed(mix);
end

function update_properties_velocities(app, mix, suffix)
    % Update velocities and Mach numbers of the mixture
    
    % Shock velocity (pre-shock or post-shock)
    switch lower(suffix)
        case {'r', '1', '_1'}
            app.(['text_u', suffix]).Value = velocity_relative(mix);
        otherwise
            app.(['text_u', suffix]).Value = mix.v_shock;
    end
    % Mach number
    app.(['text_M', suffix]).Value = app.(['text_u', suffix]).Value / soundspeed(mix);
end

function update_properties_rocket(app, mix, suffix)
    % Update rocket propellant performance parameters
    app.(['text_Aratio', suffix]).Value = mix.Aratio;
    app.(['text_Cstar', suffix]).Value = mix.cstar;
    app.(['text_Ivac', suffix]).Value = mix.I_vac;
    app.(['text_Isp', suffix]).Value = mix.I_sp;
end

function update_properties_oblique(app, mix, suffix)
    % Update rocket propellant performance parameters
    app.(['text_beta_min', suffix]).Value = mix.beta_min;
    app.(['text_beta', suffix]).Value = mix.beta;
    app.(['text_theta', suffix]).Value = mix.theta;
end

function update_error_method()
    % Update error of the method employed

end

function update_mixtures(app, results, i, FLAG_REACTANTS)
    % Update mixture tables with the results of the ith case solved
    if FLAG_REACTANTS
        species = results(i).UITable_R_Data(:, 1);
        type = results(i).UITable_R_Data(:, 4);
        temperature = results(i).UITable_R_Data(:, 5);
        ind_species = find_ind(results(i).LS, species);
        [N, Xi, ind_sort] = sort_mixture(results, 'mix1', i, ind_species);
        data = table2cell(table(species(ind_sort), Xi .* N, Xi, type(ind_sort), temperature(ind_sort)));

        app.UITable_R.Data = data;
        app.UITable_R2.Data = data(:, 1:3);
        % Update GUI: ListProducts
        app.listbox_Products.Items = results(i).LS;
    end
    species = results(i).LS';
    ind_species = find_ind(results(i).LS, species);
    [N, Xi, ind_sort] = sort_mixture(results, 'mix2', i, ind_species);
    data = table2cell(table(species(ind_sort), Xi .* N, Xi));
    app.UITable_P.Data = data;
end

function [N, Xi, ind_sort] = sort_mixture(results, mix, i, ind_species)
    % Sort given mixture in descend order
    N = results(i).(mix).N;
    [Xi, ind_sort] = sort(results(i).(mix).Xi(ind_species), 'descend');
end

function update_equivalence_ratio(app, results, i)
    % Update GUI: equivalence ratio, O/F, and percentage Fuel
    if strcmp(results(i).mix1.phi, '-')
        app.edit_phi.Value = '-';
    else
        app.edit_phi.Value = sprintf('%.5g', results(i).mix1.phi); 
    end
    app.edit_phi2.Value = app.edit_phi.Value;
    app.edit_phi3.Value = app.edit_phi.Value;
    app.edit_OF.Value = 1/results(i).mix1.FO;
    app.edit_F.Value = results(i).mix1.percentage_Fuel;
end