function [error_prop, index_max] = compute_error_prop_CEA(resultsCT, resultsCEA, varsname_x, value, varsname_y, type)
    % Compute max error of CT against CEA
    %
    % Args:
    %
    % 
    % Returns:
    %
    %
    % Examples:
    %     * 
    %     *

    % Definitions
    numCases = length(varsname_x);
    basis = 'mi';

    for i = numCases:-1:1
        varname_x = varsname_x{i};
        varname_y = varsname_y{i};

        value_CT = resultsCT.(varname_y);

        % Prepare data (units)
        switch lower(varname_y)
            case {'cp', 'cv', 'hf', 'ef', 'h', 'e', 'g', 's', 'det', 'dht', 'ds', 's0'}
                % Change units to [kJ ...]
                value_CT = value_CT * 1e-3; 
                % Property has to be divided by the basis (kg or mol)
                value_CT = value_CT ./ resultsCT.(basis);
            case {'mw', 'w'}
                % Change units to [g ...]
                value_CT = value_CT * 1e3;
        end

        try
            FLAG_CASE_CEA = resultsCEA.(varname_x) == value;
            value_CEA = resultsCEA.(varname_y)(:, FLAG_CASE_CEA);
        catch

            if strcmpi(varname_x, 'u')
                FLAG_CASE_CEA = resultsCEA.mix1.(varname_x) == value;
            else
                FLAG_CASE_CEA = resultsCEA.(type).(varname_x) == value;
            end

            value_CEA = resultsCEA.(type).(varname_y)(:, FLAG_CASE_CEA);
        end

        error_prop(i) = abs((value_CT - value_CEA) ./ value_CT);
    end
    
    % Return max error and its index from propertiesList
    [error_prop, index_max] = max(error_prop);
end