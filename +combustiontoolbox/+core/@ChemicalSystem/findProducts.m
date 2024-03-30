function [LS, ind_elements_DB] = findProducts(obj, listSpecies, varargin)
    % Find all the combinations of species from DB that can appear as
    % products for the given list of reactants
    %
    % Args:
    %     obj (Database): Database
    %     listSpecies (cell): List of reactants
    %
    % Optional Name-Value Pairs Args:
    %     * ind_elements_DB (float): Matrix NS_DB x MAX_ELEMENTS with element indeces of the species contained in the database
    %     * FLAG_BURCAT (bool): Flag indicating to look for species also in Burcat's database
    %     * FLAG_ION (bool): Flag indicating to include ionized species
    %     * FLAG_CONDENSED (bool): Flag indicating to include condensed species
    %
    % Returns:
    %     Tuple containing
    %
    %     * LS (cell): List of products
    %     * ind_elements_DB (float): Matrix NS_DB x MAX_ELEMENTS with element indeces of the species contained in the database
    %
    % Examples:
    %     * [LS, ind_elements_DB] = find_products(self, {'O2', 'N', 'eminus'})
    %     * [LS, ind_elements_DB] = find_products(self, {'O2', 'CO', 'N'}, 'flag_burcat', true)
    %     * [LS, ind_elements_DB] = find_products(self, {'O2', 'CO', 'N'}, 'flag_burcat', true, 'flag_ion', true)
    %     * [LS, ind_elements_DB] = find_products(self, {'O2', 'CO', 'N'}, 'flag_burcat', true, 'flag_ion', true, 'flag_condensed', true, 'ind', ind_elements_DB)
    %     * [LS, ind_elements_DB] = find_products(self, {'O2', 'CO', 'N'}, 'flag_burcat', true, 'flag_ion', true, 'ind', ind_elements_DB)
    
    import combustiontoolbox.utils.findIndex
    % Definitions
    MAX_ELEMENTS = 5;
    FLAG_BURCAT = obj.FLAG_BURCAT; % Flag indicating to look for species also in Burcat's database
    FLAG_ION = obj.FLAG_ION; % Flag indicating to consider ionized species
    FLAG_CONDENSED = obj.FLAG_CONDENSED; % Flag indicating to consider condensed species
    elements = combustiontoolbox.core.Elements().getElements(); % Elements list
    
    % Initialization
    FLAG_IND = false;     % Flag indicating that ind_elements_DB is an input
    LS = [];              % Initialize list of products
    LS_DB = obj.database.listSpecies; % Get list of species in the database

    % Unpack
    for i = 1:2:nargin-2
        switch lower(varargin{i})
            case {'ind', 'ind_elements', 'ind_elements_db', 'ind_db'}
                ind_elements_DB = varargin{i + 1};
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
        listSpecies{end + 1} = 'eminus';
    end

    % If FLAG_BURCAT == false, then remove the species from Third millenium database (Burcat)
    if ~FLAG_BURCAT
        LS_DB = find_species_LS(LS_DB, {}, 'any', {'_M'}, 'all');
        
        % If ind_elements_DB is given, now have to been recalculated
        if FLAG_IND
            FLAG_IND = false;
        end

    end

    % If FLAG_CONDENSED == true, then include condensed species
    if ~FLAG_CONDENSED
        NS = length(LS_DB);
        FLAG_REMOVE = false(NS, 1);
        for i = 1:NS
            if self.DB.(LS_DB{i}).phase == 1
                FLAG_REMOVE(i) = true;
            end

        end

        LS_DB(FLAG_REMOVE) = [];

        % If ind_elements_DB is given, now have to been recalculated
        if FLAG_IND
            FLAG_IND = false;
        end

    end

    % Get element indeces of the reactants
    ind_elements_R = sort(get_ind_elements(listSpecies, obj.database.species, elements, MAX_ELEMENTS), 1);

    % Get element indeces of each species in the database
    if ~FLAG_IND
        % Remove incompatible species
        LS_DB(findIndex(LS_DB, 'Air')) = [];
        % Get element indeces of each species in the database
        ind_elements_DB = get_ind_elements(LS_DB, obj.database.species, elements, MAX_ELEMENTS);
    end
    
    % Remove common vectors
    ind_elements_R = unique(ind_elements_R, 'rows');
    
    % Get cross terms 
    temp = get_cross_termms(ind_elements_R, MAX_ELEMENTS);
    
    % Join vectors
    ind_elements_R = unique([ind_elements_R; temp], 'rows');
    
    % Get list of products
    for i = length(ind_elements_R(:, 1)):-1:1
        ind_species = all(ind_elements_DB - ind_elements_R(i, :) == 0, 2);
        LS = [LS, LS_DB(ind_species)'];
    end

end

% SUB-PASS FUNCTIONS
function temp = get_cross_termms(ind_elements_R, MAX_ELEMENTS)
    % Get cross terms 

    % Initialization
    temp = [];

    % Define vector to permute
    ind_perm = sort(unique(ind_elements_R), 1, 'descend')';
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