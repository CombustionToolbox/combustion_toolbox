function app = get_current_reactants_gui(obj, app)
    % Get the species and the number of moles of the current data in
    % UITable_R

    % Get indexes of the different type of species in the mixture
    obj = gui_get_typeSpecies(obj);
    % Get current mixture
    species = obj.UITable_R.Data(:, 1);
    moles = obj.UITable_R.Data(:, 2);
    app = get_current_species_gui(obj, app, species, moles, 'S_Fuel', 'N_Fuel', 'ind_Fuel');
    app = get_current_species_gui(obj, app, species, moles, 'S_Oxidizer', 'N_Oxidizer', 'ind_Oxidizer');
    app = get_current_species_gui(obj, app, species, moles, 'S_Inert', 'N_Inert', 'ind_Inert');
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

function app = get_current_species_gui(obj, app, species, moles, species_name, moles_name, ind_name)
    % Get the species and the number of moles of the current data in
    % UITable_R for a given category (Fuel, Oxidizer, Inert)
    app.PD.(species_name) = species(obj.(ind_name))';
    try
        app.PD.(moles_name) = cell2vector(moles(obj.(ind_name)));
    catch 
        app.PD.(moles_name) = [];
    end
end