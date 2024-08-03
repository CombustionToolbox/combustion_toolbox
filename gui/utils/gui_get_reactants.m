function gui_get_reactants(app, event, varargin)
    % Get the species and the number of moles of the current data in
    % UITable_R

    % % Definitions
    % FLAG_COMPUTE_FROM_PHI = [];
    % 
    % % Unpack additional inputs
    % if nargin > 3
    %     FLAG_COMPUTE_FROM_PHI = varargin{1,1};
    % end

    % Get indexes of the different type of species in the mixture
    app = gui_get_typeSpecies(app);
    
    % Get species in the mixture
    listSpecies = app.UITable_R.Data(:, 1)';

    % Get number of moles of each species in the mixture
    moles = cell2vector(app.UITable_R.Data(:, 2))';
    
    % Define mixture
    if sum(app.indexFuel)
        set(app.mixture, listSpecies(app.indexFuel), 'fuel', moles(app.indexFuel));
    end

    if sum(app.indexOxidizer)
        set(app.mixture, listSpecies(app.indexOxidizer), 'oxidizer', moles(app.indexOxidizer));
    end

    if sum(app.indexInert)
        set(app.mixture, listSpecies(app.indexInert), 'inert', moles(app.indexInert));
    end

end

% SUB-PASS FUNCTIONS
function app = gui_get_typeSpecies(app)
    % Function that obtains the indexes of the different type of species in
    % the mixture
    typeSpecies = app.UITable_R.Data(:, 4);
    app.indexFuel     = contains(typeSpecies, 'Fuel');
    app.indexOxidizer = contains(typeSpecies, 'Oxidizer');
    app.indexInert    = contains(typeSpecies, 'Inert');
end

function self = gui_set_species(app, self, species, species_name, ind_name)
    % Get the species of the current data in
    % UITable_R for a given category (Fuel, Oxidizer, Inert)
    self.PD.(species_name) = species(app.(ind_name))';
end

function self = gui_set_species_moles_temperatures(app, self, species, moles, temperatures, species_name, moles_name, temperatures_name, ind_name)
    % Get the species and the number of moles and the temperatures of the 
    % current data in UITable_R for a given category (Fuel, Oxidizer, Inert)

    % Set species
    self = gui_set_species(app, self, species, species_name, ind_name);
    
    % Set moles
    try
        self.PD.(moles_name) = cell2vector(moles(app.(ind_name))');
    catch 
        self.PD.(moles_name) = [];
    end
    
    % Set temperatures
    try
        if ~self.Misc.FLAGS_PROP.TR
            self.PD.(temperatures_name) = cell2vector(temperatures(app.(ind_name))');
        else
            self.PD.(temperatures_name) = [];
        end
        
    catch 
        self.PD.(temperatures_name) = [];
    end

end

function moles = gui_get_moles(app, event, self, FLAG_COMPUTE_FROM_PHI)
    % Get number of moles of each species in the mixture
    if FLAG_COMPUTE_FROM_PHI
        % Set equivalence ratio
        self.PD.phi.value = gui_get_prop(self, 'phi', app.edit_phi.Value, 'direction', 'first');
        % Get proportion of oxidizers_O2
        app.PD.ratio_oxidizers_O2 = gui_ratio_oxidizers_O2(app, self.S.ind_ox_ref);
        % compute moles from equivalence ratio
        moles = gui_compute_moles_from_equivalence_ratio(app, self);
    else
        try
            if event.Indices(1, 2) ~= 3 
                moles = app.UITable_R.Data(:, 2);
            else 
                % compute moles from mole fractions
                 moles = gui_compute_moles_from_molar_fractions(app);
            end

        catch
             moles = gui_compute_moles_from_molar_fractions(app);
        end

    end

end

function moles = gui_compute_moles_from_molar_fractions(app)
    % Compute moles from mole fractions
    total_moles = sum(cell2vector(app.UITable_R.Data(:, 2)));
    mole_fractions = cell2vector(app.UITable_R.Data(:, 3));
    moles = mole_fractions * total_moles;
    if isnan(moles)
        moles = zeros(1, length(mole_fractions));
    end

end

function moles = gui_compute_moles_from_equivalence_ratio(app, self)
    % Compute moles from equivalence ratio (phi)
    
    if any(app.ind_Fuel)
        species = app.UITable_R.Data(:, 1);
        moles = app.UITable_R.Data(:, 2);
        temperatures = app.UITable_R.Data(:, 5);
        self = gui_set_species_moles_temperatures(app, self, species, moles, temperatures, 'S_Fuel', 'N_Fuel', 'T_Fuel', 'ind_Fuel');
        self = gui_set_species_moles_temperatures(app, self, species, moles, temperatures, 'S_Oxidizer', 'N_Oxidizer', 'T_Oxidizer', 'ind_Oxidizer');
        self = gui_compute_propReactants(app, self);
    end

    if any(app.ind_Oxidizer)
        moles = assign_vector2cell(moles, self.PD.phi_t/self.PD.phi.value(1) .* app.PD.ratio_oxidizers_O2, app.ind_Oxidizer);
    end

    if any(app.ind_Inert)
        moles(app.ind_Inert) = app.UITable_R.Data(app.ind_Inert, 2);
    end

end

function ratio_oxidizers_O2 = gui_ratio_oxidizers_O2(app, ind_ox_ref)
    % Get proportion of oxidizers with oxidizer of reference (default: O2)

    if ~isempty(ind_ox_ref)
        ratio_oxidizers_O2 = cell2vector(app.UITable_R.Data(app.ind_Oxidizer, 2))' / cell2vector(app.UITable_R.Data(ind_ox_ref, 2));
    else
        ratio_oxidizers_O2 = 1;
    end

end