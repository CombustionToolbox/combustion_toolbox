function ind_elements = get_ind_elements(LS, DB, elements, MAX_ELEMENTS)
    % Get element indeces of each species contained in LS
    %
    % Args:
    %     LS (cell): List of species
    %     DB (struct): Database with custom thermodynamic polynomials functions generated from NASAs 9 polynomials fits
    %     elements (cell): Elements in the periodic table
    %     MAX_ELEMENTS (float): Maximum number of elements contained in one species
    %
    % Returns:
    %     ind_elements (float): Matrix numel(LS) x MAX_ELEMENTS with element indeces of the species contained in LS
    
    % Definitions
    NS = length(LS);

    % Initialization
    ind_elements = zeros(NS, MAX_ELEMENTS);

    % Get indeces
    for i = NS:-1:1
        species = LS{i};
        temp = set_element_matrix(DB.(species).txFormula, elements);
        ind_elements(i, 1:length(temp(1, :))) = temp(1, :);
    end
    
    % Sort indeces in descend order
    ind_elements = sort(ind_elements, 2, 'descend');
end