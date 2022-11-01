function moles = convert_weight_percentage_to_moles(LS, weight_percentage, DB)
    % Convert weight percentage (wt%) to moles
    %
    % Args:
    %     LS (cell): List of species
    %     weight_percentage (float): Weight percentage of the species [%]
    %     DB (struct): Database with custom thermodynamic polynomials functions generated from NASAs 9 polynomials fits
    %
    % Returns:
    %     moles (float): Number of moles [mol]

    % Check if value is a cell
    if ~iscell(LS)
        LS = {LS};
    end

    % Get molecular weight of the species
    W = set_prop_DB(LS, 'mm', DB);
    % Convert weight percentage (wt%) to moles
    moles = weight_percentage ./ W;
end
