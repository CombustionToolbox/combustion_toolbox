function [listSpecies, indexElements_DB] = findProducts(obj, listReactants, varargin)
    % Find all the combinations of species from DB that can appear as
    % products for the given list of reactants
    %
    % Args:
    %     obj (ChemicalSystem): ChemicalSystem object
    %     listReactants (cell): List of reactants
    %
    % Optional Name-Value Pairs Args:
    %     * indexElements_DB (float): Matrix NS_DB x MAX_ELEMENTS with element indeces of the species contained in the database
    %     * FLAG_BURCAT (bool): Flag indicating to look for species also in Burcat's database
    %     * FLAG_ION (bool): Flag indicating to include ionized species
    %     * FLAG_CONDENSED (bool): Flag indicating to include condensed species
    %
    % Returns:
    %     Tuple containing
    %
    %     * listSpecies (cell): List of products
    %     * indexElements_DB (float): Matrix NS_DB x MAX_ELEMENTS with element indeces of the species contained in the database
    %
    % Examples:
    %     * [listSpecies, indexElements_DB] = findProducts(ChemicalSystem(DB), {'O2', 'N', 'eminus'})
    %     * [listSpecies, indexElements_DB] = findProducts(ChemicalSystem(DB), {'O2', 'CO', 'N'}, 'flag_burcat', true)
    %     * [listSpecies, indexElements_DB] = findProducts(ChemicalSystem(DB), {'O2', 'CO', 'N'}, 'flag_burcat', true, 'flag_ion', true)
    %     * [listSpecies, indexElements_DB] = findProducts(ChemicalSystem(DB), {'O2', 'CO', 'N'}, 'flag_burcat', true, 'flag_ion', true, 'flag_condensed', true, 'ind', indexElements_DB)
    %     * [listSpecies, indexElements_DB] = findProducts(ChemicalSystem(DB), {'O2', 'CO', 'N'}, 'flag_burcat', true, 'flag_ion', true, 'ind', indexElements_DB)
    
    % Import packages
    import combustiontoolbox.utils.findIndex

    % Definitions
    MAX_ELEMENTS = 5;
    FLAG_BURCAT = obj.FLAG_BURCAT; % Flag indicating to look for species also in Burcat's database
    FLAG_ION = obj.FLAG_ION; % Flag indicating to consider ionized species
    FLAG_CONDENSED = obj.FLAG_CONDENSED; % Flag indicating to consider condensed species
    elements = combustiontoolbox.core.Elements().getElements(); % Elements list
    
    % Initialization
    FLAG_IND = false; % Flag indicating that indexElements_DB is an input
    listSpecies = []; % Initialize list of products
    listSpecies_DB = obj.database.listSpecies; % Get list of species in the database

    % Unpack
    for i = 1:2:nargin-2
        switch lower(varargin{i})
            case {'ind', 'ind_elements', 'ind_elements_db', 'ind_db', 'indexelements_db'}
                indexElements_DB = varargin{i + 1};
                FLAG_IND = true;
            case {'flag', 'flag_burcat'}
                FLAG_BURCAT = varargin{i + 1};
            case {'flag_ion', 'flag_ions'}
                FLAG_ION = varargin{i + 1};
            case {'flag_condensed'}
                FLAG_CONDENSED = varargin{i + 1};
        end
        
    end
    
    % If FLAG_ION == true, then find ionized species
    if FLAG_ION
        listReactants{end + 1} = 'eminus';
    end

    % If FLAG_BURCAT == false, then remove the species from Third millenium database (Burcat)
    if ~FLAG_BURCAT
        listSpecies_DB = findSpecies(listSpecies_DB, {}, 'any', {'_M'}, 'all');
        
        % If indexElements_DB is given, now have to been recalculated
        if FLAG_IND
            FLAG_IND = false;
        end

    end

    % If FLAG_CONDENSED == true, then include condensed species
    if ~FLAG_CONDENSED
        NS = length(listSpecies_DB);
        FLAG_REMOVE = false(NS, 1);
        for i = 1:NS
            if self.DB.(listSpecies_DB{i}).phase == 1
                FLAG_REMOVE(i) = true;
            end

        end

        listSpecies_DB(FLAG_REMOVE) = [];

        % If indexElements_DB is given, now have to been recalculated
        if FLAG_IND
            FLAG_IND = false;
        end

    end

    % Get element indeces of the reactants
    indexElements = sort(getIndexElements(listReactants, obj.database.species, elements, MAX_ELEMENTS), 1);

    % Get element indeces of each species in the database
    if ~FLAG_IND
        % Remove incompatible species
        listSpecies_DB(findIndex(listSpecies_DB, 'Air')) = [];
        % Get element indeces of each species in the database
        indexElements_DB = getIndexElements(listSpecies_DB, obj.database.species, elements, MAX_ELEMENTS);
    end
    
    % Remove common vectors
    indexElements = unique(indexElements, 'rows');
    
    % Get cross terms 
    temp = getCrossTermms(indexElements, MAX_ELEMENTS);
    
    % Join vectors
    indexElements = unique([indexElements; temp], 'rows');
    
    % Get list of products
    for i = length(indexElements(:, 1)):-1:1
        indexSpecies = all(indexElements_DB - indexElements(i, :) == 0, 2);
        listSpecies = [listSpecies, listSpecies_DB(indexSpecies)'];
    end

end

% SUB-PASS FUNCTIONS
function temp = getCrossTermms(indexElements, MAX_ELEMENTS)
    % Get cross terms 

    % Initialization
    temp = [];

    % Define vector to permute
    ind_perm = sort(unique(indexElements), 1, 'descend')';
    try
        ind_perm = [ind_perm, zeros(1, MAX_ELEMENTS - 2)];
    catch
        ind_perm = [ind_perm', zeros(1, MAX_ELEMENTS - 2)];
    end

    % Get length
    N_elements = length(ind_perm);

    % Get permutations [element1, element2, element3, element4, element5]
    for i = 1:N_elements
        if ind_perm(i) == 0, continue; end
        
        for j = i + 1:N_elements

            for k = j + 1:N_elements

                for l = k + 1:N_elements-1
                    temp_add = [ind_perm(i), ind_perm(j), ind_perm(k), ind_perm(l), 0];
                    temp = [temp; temp_add];
                end

            end

        end

    end

end