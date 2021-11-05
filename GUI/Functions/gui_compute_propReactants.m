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
    if ~isempty(app.PD.S_Fuel) && ~isempty(app.PD.S_Oxidizer)
        % Compute percentage Fuel, Oxidant/Fuel ratio and equivalence ratio
        R = app.PD.R_Fuel + app.PD.R_Oxidizer + app.PD.R_Inert;
        app.PS.strR{1}.percentage_Fuel = sum(app.PD.R_Fuel(:, 1)) / sum(R(:, 1)) * 100;
        app.PS.strR{1}.FO = sum(app.PD.R_Fuel(:, 1)) / sum(app.PD.R_Oxidizer(:, 1));
        app.PS.strR{1}.FO_st = sum(app.PD.R_Fuel(:, 1)) / (app.PS.strR_Fuel.x + app.PS.strR_Fuel.y/4 - app.PS.strR_Fuel.z/2);
        app.PS.strR{1}.OF = 1/app.PS.strR{1}.FO;
        app.PS.strR{1}.phi = app.PS.strR{1}.FO / app.PS.strR{1}.FO_st;
    else
        app.PS.strR{1}.percentage_Fuel = 0;
        app.PS.strR{1}.FO = inf;
        app.PS.strR{1}.OF = 0;
        app.PS.strR{1}.phi = '-';
    end
end