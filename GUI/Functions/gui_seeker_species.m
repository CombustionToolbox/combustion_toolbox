function species = gui_seeker_species(obj, event)
    % Return the species that match with the species introduced in the
    % finder.
    try
        seekSpecies = event.Value;
        listSpecies = obj.S.LS_DB;
        index = 0;
        seekIndex = false;
        while ~seekIndex
            index = index + 1;
            seekIndex = startsWith(listSpecies{index}, seekSpecies, 'IgnoreCase', false);
        end
        if index
            species = listSpecies{index}; % Species found in the Database
        else
            species = []; % Species not found in the Database
        end
    catch ME
      errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
      ME.stack(1).name, ME.stack(1).line, ME.message);
      fprintf('%s\n', errorMessage);
      uiwait(warndlg(errorMessage));
    end
end