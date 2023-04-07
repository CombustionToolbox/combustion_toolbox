function self = gui_compute_propReactants(app, self)
    % Function that compute fundamental properties (e.g., equivalence ratio,
    % molar fractions, ...) of the mixture
    %
    % Args:
    %     app (object): Combustion Toolbox app object
    %     self (struct): Data of the mixture, conditions, and databases
    %
    % Returns:
    %     self (struct): Data of the mixture, conditions, and databases
    
    % Get temperature and pressure of the mixture
    self.PD.TR.value = gui_get_prop(app, 'TR', app.PR1.Value);
    self.PD.pR.value = gui_get_prop(app, 'pR', app.PR2.Value);

    % Flag if the number of moles of fuel, oxidant and inert species
    % is specified. If not, consider 1 mole for the fuel and calculate
    % the remaining moles from the equivalence relation
    self = get_FLAG_N(self);

    % Compute properties for the first temperature value (in case PR1 is a vector)
    self.PD.TR.value = self.PD.TR.value(1);
    
    % Compute stoichiometric matrix and properties of the mixture
    self = define_FOI(self, 1);
end