function gui_write_results(obj, results, i, varargin)
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
    update_properties(obj, results, i);
    % Update mixture (tables)
    update_mixtures(obj, results, i, FLAG_REACTANTS);
    % Update equivalence ratio
    update_equivalence_ratio(obj, results, i);
end

% SUB-PASS FUNCTIONS
function update_properties(obj, results, i)
    mix1 = results(i).mix1;
    mix2 = results(i).mix2;
    obj.text_TR.Value = temperature(mix1);
    obj.text_TP.Value = temperature(mix2);
    obj.text_pR.Value = pressure(mix1);
    obj.text_pP.Value = pressure(mix2);
    obj.text_rR.Value = density(mix1);
    obj.text_rP.Value = density(mix2);
    obj.text_hR.Value = enthalpy_mass(mix1);
    obj.text_hP.Value = enthalpy_mass(mix2);
    obj.text_eR.Value = intEnergy_mass(mix1);
    obj.text_eP.Value = intEnergy_mass(mix2);
    obj.text_cpR.Value = cp_mass(mix1);
    obj.text_cpP.Value = cp_mass(mix2);
    try
        obj.text_sR.Value = entropy_mass(mix1);
    catch
        obj.text_sR.Value = Inf;
    end
    obj.text_sP.Value = entropy_mass(mix2);
    obj.text_gammaR.Value = adiabaticIndex(mix1);
    obj.text_gammaP.Value = adiabaticIndex(mix2);
    obj.text_WR.Value = MolecularWeight(mix1);
    obj.text_WP.Value = MolecularWeight(mix2);
    obj.text_soundR.Value = soundspeed(mix1);
    obj.text_soundP.Value = soundspeed(mix2);
    if strcmpi(results(i).ProblemType, 'TP')
        obj.text_error_problem.Value = mix2.error_moles;
    else
        obj.text_error_problem.Value = mix2.error_problem;
    end
    if contains(results(i).ProblemType, 'SHOCK', 'IgnoreCase', true) || contains(results(i).ProblemType, 'DET', 'IgnoreCase', true) || contains(results(i).ProblemType, 'ROCKET', 'IgnoreCase', true)
        obj.text_uR.Value = velocity_relative(mix1);
        obj.text_uP.Value = mix2.v_shock;
        obj.text_MR.Value = velocity_relative(mix1)/soundspeed(mix1);
        obj.text_MP.Value = mix2.v_shock/soundspeed(mix2);
        if contains(results(i).ProblemType, 'ROCKET', 'IgnoreCase', true)
            obj.text_Aratio_2.Value = mix2.Aratio;
            
            if isfield(results(i), 'mix2_c')
                mix3 = results(i).mix2_c;
                FLAG_IAC = false;
            else
                mix3 = results(i).mix3;
                FLAG_IAC = true;
            end
            obj.text_TP_3.Value = temperature(mix3);
            obj.text_pP_3.Value = pressure(mix3);
            obj.text_rP_3.Value = density(mix3);
            obj.text_hP_3.Value = enthalpy_mass(mix3);
            obj.text_eP_3.Value = intEnergy_mass(mix3);
            obj.text_cpP_3.Value = cp_mass(mix3);
            obj.text_sP_3.Value = entropy_mass(mix3);
            obj.text_gammaP_3.Value = adiabaticIndex(mix3);
            obj.text_WP_3.Value = MolecularWeight(mix3);
            obj.text_soundP_3.Value = soundspeed(mix3);
            obj.text_uP_3.Value = mix3.v_shock;
            obj.text_MP_3.Value = mix3.v_shock/soundspeed(mix3);
            obj.text_Aratio_3.Value = mix3.Aratio;
            obj.text_Cstar_3.Value = mix3.cstar;
            obj.text_Ivac_3.Value = mix3.I_vac;
            obj.text_Isp_3.Value = mix3.I_sp;
            
            if ~isempty(results(i).strP) || ~FLAG_IAC

                if ~FLAG_IAC
                    mix4 = results(i).mix3;
                else
                    mix4 = results(i).strP;
                end
                
                obj.text_TP_4.Value = temperature(mix4);
                obj.text_pP_4.Value = pressure(mix4);
                obj.text_rP_4.Value = density(mix4);
                obj.text_hP_4.Value = enthalpy_mass(mix4);
                obj.text_eP_4.Value = intEnergy_mass(mix4);
                obj.text_cpP_4.Value = cp_mass(mix4);
                obj.text_sP_4.Value = entropy_mass(mix4);
                obj.text_gammaP_4.Value = adiabaticIndex(mix4);
                obj.text_WP_4.Value = MolecularWeight(mix4);
                obj.text_soundP_4.Value = soundspeed(mix4);
                obj.text_uP_4.Value = mix4.v_shock;
                obj.text_MP_4.Value = mix4.v_shock/soundspeed(mix4);
                obj.text_Aratio_4.Value = mix4.Aratio;
                obj.text_Cstar_4.Value = mix4.cstar;
                obj.text_Ivac_4.Value = mix4.I_vac;
                obj.text_Isp_4.Value = mix4.I_sp;

                % Get max error method
                obj.text_error_problem.Value = max([mix2.error_problem, mix3.error_problem, mix4.error_problem]);
            else
                obj.text_TP_4.Value = 0;
                obj.text_pP_4.Value = 0;
                obj.text_rP_4.Value = 0;
                obj.text_hP_4.Value = 0;
                obj.text_eP_4.Value = 0;
                obj.text_cpP_4.Value = 0;
                obj.text_sP_4.Value = 0;
                obj.text_gammaP_4.Value = 0;
                obj.text_WP_4.Value = 0;
                obj.text_soundP_4.Value = 0;
                obj.text_uP_4.Value = 0;
                obj.text_MP_4.Value = 0;
                obj.text_Aratio_4.Value = 0;
                obj.text_Cstar_4.Value = 0;
                obj.text_Ivac_4.Value = 0;
                obj.text_Isp_4.Value = 0;

                % Get max error method
                obj.text_error_problem.Value = max(mix2.error_problem, mix3.error_problem);
            end
        end
    end
end

function update_mixtures(obj, results, i, FLAG_REACTANTS)
    % Update mixture tables with the results of the ith case solved
    if FLAG_REACTANTS
        species = results(i).UITable_R_Data(:, 1);
        type = results(i).UITable_R_Data(:, 4);
        temperature = results(i).UITable_R_Data(:, 5);
        ind_species = find_ind(results(i).LS, species);
        [N, Xi, ind_sort] = sort_mixture(results, 'mix1', i, ind_species);
        data = table2cell(table(species(ind_sort), Xi .* N, Xi, type, temperature));

        obj.UITable_R.Data = data;
        obj.UITable_R2.Data = data(:, 1:3);
        % Update GUI: ListProducts
        obj.listbox_Products.Items = results(i).LS;
    end
    species = results(i).LS';
    ind_species = find_ind(results(i).LS, species);
    [N, Xi, ind_sort] = sort_mixture(results, 'mix2', i, ind_species);
    data = table2cell(table(species(ind_sort), Xi .* N, Xi));
    obj.UITable_P.Data = data;
end

function [N, Xi, ind_sort] = sort_mixture(results, mix, i, ind_species)
    % Sort given mixture in descend order
    N = results(i).(mix).N;
    [Xi, ind_sort] = sort(results(i).(mix).Xi(ind_species), 'descend');
end

function update_equivalence_ratio(obj, results, i)
    % Update GUI: equivalence ratio, O/F, and percentage Fuel
    if strcmp(results(i).mix1.phi, '-')
        obj.edit_phi.Value = '-';
    else
        obj.edit_phi.Value = sprintf('%.5g', results(i).mix1.phi); 
    end
    obj.edit_phi2.Value = obj.edit_phi.Value;
    obj.edit_phi3.Value = obj.edit_phi.Value;
    obj.edit_OF.Value = 1/results(i).mix1.FO;
    obj.edit_F.Value = results(i).mix1.percentage_Fuel;
end