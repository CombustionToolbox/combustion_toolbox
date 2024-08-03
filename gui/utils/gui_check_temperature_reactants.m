function [app, temperature, FLAG_FIXED] = gui_check_temperature_reactants(app, listSpecies, temperature, numSpecies)
    % Check if is a condensed species with a fixed temperature
    FLAG_FIXED = false;
    
    
    for i = 1:numSpecies
        % Get Species object
        species = app.database.species.(listSpecies{i});
        
        % Check if condensed species can be only evaluated a particular temperature
        if isscalar(species.T)
            temperature{i} = app.database.species.(listSpecies{i}).T;
            FLAG_FIXED = true;
        end

    end

    if FLAG_FIXED
        app.PR1.Value = sprintf('%.4g', max(cell2mat(temperature))); 
    end

end