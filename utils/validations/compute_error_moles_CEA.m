function error_moles = compute_error_moles_CEA(resultsCT, resultsCEA, varname_x, value, varname_y, listSpecies, tolMoles)
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
    
    % Import packages
    import combustiontoolbox.utils.findIndex
    
    % Get index species
    index_species_CT = findIndex(resultsCT.chemicalSystem.listSpecies, listSpecies);
    index_case_CEA = find(round(resultsCEA.(varname_x), 4) == value);

    moles_CT = resultsCT.(varname_y)(index_species_CT);
    moles_CEA = resultsCEA.(varname_y)(:, index_case_CEA);
    
    % 
    moles_CEA = moles_CEA(:);

    % if any(size(moles_CT) ~= size(moles_CEA))
    %     moles_CEA = moles_CEA';
    % end

    index_check = find(moles_CT > 10^floor(log10(tolMoles) + 2));

    [error_moles, index] = max(abs((moles_CT(index_check) - moles_CEA(index_check)) ./ moles_CT(index_check)));
end