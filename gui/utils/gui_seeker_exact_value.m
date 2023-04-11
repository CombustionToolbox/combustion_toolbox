function value = gui_seeker_exact_value(app, event, ListValues)
    % Return the value that match with the value introduced in the finder
    try
        seekValue = event.Value;
        index = 0;
        seekIndex = false;
        while ~seekIndex
            index = index + 1;
            seekIndex = strcmp(ListValues{index}, seekValue);
        end

        if index
            value = ListValues{index}; % Value found in ListValues
        else
            value = []; % Value not found in ListValues
        end

    catch ME
        errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
        ME.stack(1).name, ME.stack(1).line, ME.message);
        fprintf('%s\n', errorMessage);
    end
    
end