function indexElements = getIndexElements(listSpecies, database, elements, MAX_ELEMENTS)
    % Get element indeces of each species contained in LS
    %
    % Args:
    %     listSpecies (cell): List of species
    %     database (Database): Database with custom thermodynamic polynomials functions generated from NASAs 9 polynomials fits
    %     elements (cell): Elements in the periodic table
    %     MAX_ELEMENTS (float): Maximum number of elements contained in one species
    %
    % Returns:
    %     indexElements (float): Matrix numel(LS) x MAX_ELEMENTS with element indeces of the species contained in LS
    %
    % Example:
    %     indexElements = getIndexElements(listSpecies, database, elements, MAX_ELEMENTS)
    
    % Definitions
    NS = length(listSpecies);

    % Initialization
    indexElements = zeros(NS, MAX_ELEMENTS);

    % Get indeces
    for i = NS:-1:1
        species = listSpecies{i};
        temp = database.(species).getElementMatrix(elements);
        indexElements(i, 1:length(temp(1, :))) = temp(1, :);
    end
    
    % Sort indeces in descend order
    indexElements = sort(indexElements, 2, 'descend');
end