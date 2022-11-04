function [app, temperature, FLAG_FIXED] = gui_check_temperature_reactants(app, DB, species, temperature, Nspecies)
    % Check if is a condensed species with a fixed temperature
    FLAG_FIXED = true;
    for i = 1:Nspecies
        if length(DB.(species{i}).T) == 1
            temperature{i} = DB.(species{i}).T;
        else
            FLAG_FIXED = false;
        end
    end

    if FLAG_FIXED
        app.PR1.Value = sprintf('%.4g', max(cell2mat(temperature))); 
    end
end