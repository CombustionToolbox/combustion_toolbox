function self = gui_create_temp_app(app, event, FLAG_COMPUTE_FROM_PHI)
    % Function that creates an self required for preliminary calculations
    
    % Initialize self (fast: transfer DB)
    self = App('fast', app.DB_master, app.DB);
    % Get reactant species
    self = gui_get_reactants(app, event, self, FLAG_COMPUTE_FROM_PHI);
    % Compute properties of the mixture
    self = gui_compute_propReactants(app, self);
end