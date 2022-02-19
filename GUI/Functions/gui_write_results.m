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
    obj.text_sR.Value = entropy_mass(mix1);
    obj.text_sP.Value = entropy_mass(mix2);
    obj.text_gammaR.Value = adiabaticIndex(mix1);
    obj.text_gammaP.Value = adiabaticIndex(mix2);
    obj.text_WR.Value = meanMolecularWeight(mix1);
    obj.text_WP.Value = meanMolecularWeight(mix2);
    obj.text_soundR.Value = soundspeed(mix1);
    obj.text_soundP.Value = soundspeed(mix2);
    obj.text_q.Value = obj.text_hP.Value - obj.text_hR.Value;
    obj.text_error_moles.Value = mix2.error_moles;
    if strcmpi(results(i).ProblemType, 'TP')
        obj.text_error_problem.Value = mix2.error_moles;
    else
        obj.text_error_problem.Value = mix2.error_problem;
    end
    if contains(results(i).ProblemType, 'SHOCK', 'IgnoreCase', true) || contains(results(i).ProblemType, 'DET', 'IgnoreCase', true)
        obj.text_uR.Value = velocity_relative(mix1);
        obj.text_uP.Value = mix2.v_shock;
        obj.text_MR.Value = velocity_relative(mix1)/soundspeed(mix1);
        obj.text_MP.Value = mix2.v_shock/soundspeed(mix2);
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
    obj.edit_OF.Value = 1/results(i).mix1.FO;
    obj.edit_F.Value = results(i).mix1.percentage_Fuel;
end