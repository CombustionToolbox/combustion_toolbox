function app = gui_create_temp_app(obj, event, FLAG_COMPUTE_FROM_PHI)
    % Function that creates an app required for preliminary calculations
    
    % Initialize app (fast: transfer DB)
    app = App('fast', obj.DB_master, obj.DB);
    % Get reactant species
    app = gui_get_reactants(obj, event, app, FLAG_COMPUTE_FROM_PHI);
    % Compute properties of the mixture
    app = gui_compute_propReactants(obj, app);
end