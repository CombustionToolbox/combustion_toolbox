function error_moles = compute_error_moles_CEA(results1, results2, varname_x, value, varname_y, species)
    % Compute max error of CT against CEA
    dataname_y = get_dataname(varname_y);
    results1.(varname_y) = cell2vector(select_data(results1, dataname_y), varname_y);
    index_species_CT = find_ind(results1.Misc.LS_original, species);
    index_case_CEA = results2.(varname_x) == value;
    
    moles_CT = results1.(varname_y)(index_species_CT);
    moles_CEA = results2.(varname_y)(:, index_case_CEA);
    index_check = moles_CT > results1.TN.tolN;
    
    error_moles = max(abs((moles_CT(index_check) - moles_CEA(index_check)) ./ moles_CT(index_check)));
end

% SUB-PASS FUNCTIONS
function dataname = get_dataname(var)
    if strcmpi(var, 'phi')
        dataname = 'PD.phi.value';
    else
        dataname = 'PS.strP';
    end
end

function dataselected = select_data(self, dataname)
    index = strfind(dataname, '.');
    index = [index, length(dataname) + 1];
    N = length(index);
    dataselected = self;
    pos1 = 1;
    for i = 1:N
        pos2 = index(i) - 1;
        varname = dataname(pos1:pos2);
        dataselected = dataselected.(varname);
        pos1 = index(i) + 1;
    end
end