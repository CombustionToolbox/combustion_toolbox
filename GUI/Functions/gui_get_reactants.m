function app = gui_get_reactants(obj, event, app)
    % Get the species and the number of moles of the current data in
    % UITable_R

    % Get indexes of the different type of species in the mixture
    obj = gui_get_typeSpecies(obj);
    % Get species in the mixture
    species = obj.UITable_R.Data(:, 1);
    % Get number of moles of each species in the mixture
    moles = gui_get_moles(obj, event);
    % Set mixture into variable app
    app = gui_set_species_moles(obj, app, species, moles, 'S_Fuel', 'N_Fuel', 'ind_Fuel');
    app = gui_set_species_moles(obj, app, species, moles, 'S_Oxidizer', 'N_Oxidizer', 'ind_Oxidizer');
    app = gui_set_species_moles(obj, app, species, moles, 'S_Inert', 'N_Inert', 'ind_Inert');
    % Compute proportion inerts/O2
    app = compute_proportion_inerts_O2(app);
end

% SUB-PASS FUNCTIONS
function obj = gui_get_typeSpecies(obj)
    % Function that obtains the indexes of the different type of species in
    % the mixture
    typeSpecies = obj.UITable_R.Data(:, 4);
    obj.ind_Fuel     = contains(typeSpecies, 'Fuel');
    obj.ind_Oxidizer = contains(typeSpecies, 'Oxidizer');
    obj.ind_Inert    = contains(typeSpecies, 'Inert');
end

function app = gui_set_species_moles(obj, app, species, moles, species_name, moles_name, ind_name)
    % Get the species and the number of moles of the current data in
    % UITable_R for a given category (Fuel, Oxidizer, Inert)
    app.PD.(species_name) = species(obj.(ind_name))';
    try
        app.PD.(moles_name) = cell2vector(moles(obj.(ind_name)));
    catch 
        app.PD.(moles_name) = [];
    end
end

function moles = gui_get_moles(obj, event)
    % Get number of moles of each species in the mixture
    if event.Indices(1, 2) ~= 3 
        moles = obj.UITable_R.Data(:, 2);
    else % compute moles from mole fractions
        total_moles = sum(cell2vector(obj.UITable_R.Data(:, 2)));
        mole_fractions = cell2vector(obj.UITable_R.Data(:, 3));
        moles = mole_fractions * total_moles;
    end
end