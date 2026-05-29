function [error_prop, index_max] = compute_error_prop_Cantera(resultsCT, resultsCantera, varsname_x, value, varsname_y, type)
    % Compute max error of CT against Cantera
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
            FLAG_CASE_Cantera = abs(resultsCantera.(varname_x) - value) < 1e-8;
            value_Cantera = resultsCantera.(varname_y)(FLAG_CASE_Cantera);
        catch

            if strcmpi(varname_x, 'u')
                FLAG_CASE_Cantera = resultsCantera.mix1.(varname_x) == value;
            else
                FLAG_CASE_Cantera = abs(resultsCantera.(varname_x) - value) < 1e-8;
            end

            value_Cantera = resultsCantera.(type).(varname_y)(FLAG_CASE_Cantera);
        end

        error_prop(i) = abs((value_CT - value_Cantera) ./ value_CT);
    end
    
    % Return max error and its index from propertiesList
    [error_prop, index_max] = max(error_prop);
end