function error_moles = compute_error_prop_CEA(results1, results2, varsname_x, value, varsname_y, type)
    % Compute max error of CT against CEA
    for i = length(varsname_x):-1:1
        varname_x = varsname_x{i};
        varname_y = varsname_y{i};

        dataname_y = get_dataname(varname_y, type);
        value_CT = cell2vector(select_data(results1, dataname_y), varname_y);
        try
            index_case_CEA = results2.(varname_x) == value;
            value_CEA = results2.(varname_y)(:, index_case_CEA);
        catch
            if strcmpi(varname_x, 'u')
                index_case_CEA = results2.mix1.(varname_x) == value;
            else
                index_case_CEA = results2.(type).(varname_x) == value;
            end
            value_CEA = results2.(type).(varname_y)(:, index_case_CEA);
        end

        error_moles(i) = abs((value_CT - value_CEA) ./ value_CT);
    end
    error_moles = max(error_moles);
end

% SUB-PASS FUNCTIONS
function dataname = get_dataname(var, type)
    if strcmpi(var, 'phi')
        dataname = 'PD.phi.value';
    elseif strcmpi(type, 'mix1')
        dataname = 'PS.strR';
    elseif strcmpi(type, 'mix2')
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