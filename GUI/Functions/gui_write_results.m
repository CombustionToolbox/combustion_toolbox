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
end

% SUB-PASS FUNCTIONS
function update_properties(obj, results, i)
    mix1 = results.mix1{i};
    mix2 = results.mix2{i};
    obj.text_TR.Value = temperature(mix1);
    obj.text_TP.Value = temperature(mix2);
    obj.text_pR.Value = pressure(mix1);
    obj.text_pP.Value = pressure(mix1);
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
    if sscanf(results.ProblemType, '%f') > 6
        obj.text_uR.Value = velocity_relative(mix1);
        obj.text_uP.Value = velocity_relative(mix2);
        obj.text_MR.Value = velocity_relative(mix1)/soundspeed(mix1);
        obj.text_MP.Value = velocity_relative(mix2)/soundspeed(mix2);
    end
end

function update_mixtures(obj, results, i, FLAG_REACTANTS)
    % Update mixture tables with the results of the ith case solved
    if FLAG_REACTANTS
        species = results.UITable_R_Data(:, 1);
        type = results.UITable_R_Data(:, 4);
        temperature = results.UITable_R_Data(:, 5);
        ind_species = find_ind(results.LS_products, species);
        [N, Xi, ind_sort] = sort_mixture(results, 'mix1', i, ind_species);
        data = table2cell(table(species(ind_sort), Xi .* N, Xi, type, temperature));
        obj.UITable_R.Data = data;
        obj.UITable_R2.Data = data(:, 1:3);
    end
    species = results.LS';
    ind_species = find_ind(results.LS, species);
    [N, Xi, ind_sort] = sort_mixture(results, 'mix2', i, ind_species);
    data = table2cell(table(species(ind_sort), Xi .* N, Xi));
    obj.UITable_P.Data = data;
end

function [N, Xi, ind_sort] = sort_mixture(results, mix, i, ind_species)
    % Sort given mixture in descend order
    N = results.(mix){i}.N;
    [Xi, ind_sort] = sort(results.(mix){i}.Xi(ind_species), 'descend');
end