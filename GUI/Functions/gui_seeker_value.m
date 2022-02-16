function value = gui_seeker_value(obj, event, ListValues)
    % Return the value that match with the value introduced in the finder
    try
        seekValue = event.Value;
        if isempty(seekValue)
            value = [];
            return
        end
        for i = length(ListValues):-1:1
            seekIndex(i) = startsWith(ListValues{i}, seekValue, 'IgnoreCase', false);
        end
        if any(seekIndex)
            value = ListValues(seekIndex); % Value found in ListValues
        else
            value = []; % Value not found in ListValues
        end
    catch ME
      errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
      ME.stack(1).name, ME.stack(1).line, ME.message);
      fprintf('%s\n', errorMessage);
      uiwait(warndlg(errorMessage));
    end
end