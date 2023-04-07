function self = gui_create_temp_app(app, event, FLAG_COMPUTE_FROM_PHI)
    % Function that creates a self struct required for preliminary calculations
    % 
    % Args:
    %     app (object): Combustion Toolbox app object
    %     event (object): Event object
    %     FLAG_COMPUTE_FROM_PHI (logical): Flag to compute properties from the equivalence ratio
    %
    % Returns:
    %     self (struct): Struct containing the properties of the mixture and the databases

    % Initialize self (fast: transfer DB)
    self = App('fast', app.DB_master, app.DB);
    % Get reactant species
    self = gui_get_reactants(app, event, self, FLAG_COMPUTE_FROM_PHI);
    % Compute properties of the mixture
    self = gui_compute_propReactants(app, self);
end