function LS = find_products(self, species, varargin)
    % Find all the combinations of species from DB that can appear as
    % products for the given list of reactants
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %     species (cell): List of reactants
    %
    % Optional Args:
    %     flag (bool): Flag to consider Burcat's database (Third millenium)
    %
    % Returns:
    %     LS (cell): List of products
    %
    % Examples:
    %     LS = find_products({'O2'}, DB)
    %     LS = find_products({'O2', 'CO', 'N'}, DB)
    %     LS = find_products({'O2', 'N', 'eminus'}, DB)
    %     LS = find_products({'O2', 'CO', 'N'}, DB, 'Flag', true)
    
    % Definitions
    MAX_ELEMENTS = 5;
    [elements, ~] = set_elements(); % Elements list
    
    % Initialization
    FLAG_BURCAT = false;
    LS = [];
    % Unpack
    for i = 1:2:nargin-2
        switch lower(varargin{i})
            case {'flag', 'flag_burcat'}
                FLAG_BURCAT = varargin{i + 1};
        end
        
    end

    % If FLAG_BURCAT == true, then remove the species from Third millenium database
    if ~FLAG_BURCAT
        self.S.LS_DB = find_species_LS(self.S.LS_DB, {}, 'any', {'_M'}, 'all');
    end

    % Remove incompatible species
    self.S.LS_DB(find_ind(self.S.LS_DB,'Air')) = [];
    
    % Get element indeces of the reactants
    ind_elements_R = sort(sort(get_elements(species, self.DB, elements, MAX_ELEMENTS), 2, 'descend'), 1);

    % Get element indeces of each species in the database
    ind_elements_DB = sort(get_elements(self.S.LS_DB, self.DB, elements, MAX_ELEMENTS), 2, 'descend');
    
    % Remove common vectors
    ind_elements_R = unique(ind_elements_R, 'rows');
    
    % Get cross terms 
    temp = get_cross_termms(ind_elements_R, MAX_ELEMENTS);
    
    % Join vectors
    ind_elements_R = unique([ind_elements_R; temp], 'rows');
    
    % Get list of products
    for i = length(ind_elements_R(:, 1)):-1:1
        ind_species = all(ind_elements_DB - ind_elements_R(i, :) == 0, 2);
        LS = [LS, self.S.LS_DB(ind_species)'];
    end

end

% SUB-PASS FUNCTIONS
function ind_elements = get_elements(LS, DB, elements, MAX_ELEMENTS)
    % Get element indeces of each species contained in LS
    
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
    
end

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