function self = get_oxidizer_reference(self, varargin)
    % Get oxidizer of reference for computations with the equivalence ratio
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %
    % Returns:
    %     self (struct): Data of the mixture, conditions, and databases, 
    %                   included the oxidizer of reference which can be
    %                   obtained as self.S.ind_ox_ref
    
    % Definitions
    LS = self.S.LS;

    % Unpack
    if nargin > 1
        LS = {}; self.PD.S_Oxidizer = {};
        self.PD.S_Oxidizer(1, :) = varargin{1};
        LS(1, :) = varargin{1};
    end

    % Check if there are oxidizers in the mixtures
    if isempty(self.PD.S_Oxidizer)
        self.S.ind_ox_ref = [];
        self.S.atoms_ox_ref = NaN;
        return
    end

    % If O2 or O2(L) are included as oxidizers these species will be
    % selected as reference oxidizers in this order. Otherwise, the first
    % oxidizer with oxygen as element will be selected.
    if any(ismember(self.PD.S_Oxidizer, 'O2'))
        self.S.ind_ox_ref = find_ind(LS, 'O2');
        self.S.atoms_ox_ref = 2;
    elseif any(ismember(self.PD.S_Oxidizer, 'O2bLb'))
        self.S.ind_ox_ref = find_ind(LS, 'O2bLb');
        self.S.atoms_ox_ref = 2;
    else
        % Get first oxidizer with oxygen as element
        temp_ind = find(contains(self.PD.S_Oxidizer, 'O'), 1);
        species = self.PD.S_Oxidizer{temp_ind};
        % Find index of reference oxidizer
        self.S.ind_ox_ref = find_ind(LS, species);
        % Find position oxygen element
        temp_ind_O = find(self.PD.S_Oxidizer{temp_ind} == 'O');
        % Get position numbers and letters
        [temp_ind_1, temp_ind_2] = regexp(species, '\w\d*');
        % Find position oxygen element in the temp variable index
        temp_ind = find(temp_ind_1 == temp_ind_O);
        % Set number of elements of oxygen in the reference oxidizer
        self.S.atoms_ox_ref = sscanf(species(temp_ind_1(temp_ind) + 1:temp_ind_2(temp_ind)), '%f');
    end

end