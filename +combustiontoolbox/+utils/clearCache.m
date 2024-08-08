function clearCache()
    % Function to clear the variables stored in the cache
    %
    % Example:
    %     clearCache()

    % Definitions
    listFunctions = getListFunctions();

    % Clear cache
    for i = 1:length(listFunctions)
        clear(listFunctions{i});
    end

end

% SUB-PASS FUNCTIONS
function listFunctions = getListFunctions()
    % Function to get the list of functions to clear
    %
    % Example:
    %     getListFunctions()

    listFunctions = {
        'species_cP',...
        'species_DeT',...
        'species_DhT',...
        'species_g0',...
        'species_h0',...
        'species_s0',...
        'Database',...
        'NasaDatabase'
    };
end