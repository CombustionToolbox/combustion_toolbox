function app = gui_compute_propReactants(obj, app)
    % Function that compute fundamental properties (e.g., equivalence ratio,
    % molar fractions, ...) of the mixture.
    
    % Get temperature and pressure of the mixture
    app.PD.TR.value = gui_get_prop(obj, 'TR', obj.PR1.Value);
    app.PD.pR.value = gui_get_prop(obj, 'pR', obj.PR2.Value);
    % Flag if the number of moles of fuel, oxidant and inert species
    % is specified. If not, consider 1 mole for the fuel and calculate
    % the remaining moles from the equivalence relation.
    app = get_FLAG_N(app);
    % Compute stoichiometric matrix and properties of the mixture
    app = Define_FOI(app, 1);
end