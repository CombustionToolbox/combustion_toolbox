function gui_get_reactants(app)
    % FUEL
    if sum(app.ind_Fuel)
        app.PD.S_Fuel = regexprep(app.UITable_R.Data(app.ind_Fuel,1),'\s+','')';
        app.PD.N_Fuel = app.UITable_R.Data{app.ind_Fuel, 2};
    else
        app.PD.S_Fuel = [];
    end
    % OXIDIZER
    if sum(app.ind_Oxidizer)
        app.PD.S_Oxidizer = regexprep(app.UITable_R.Data(app.ind_Oxidizer,1),'\s+','')';
        app.PD.N_Oxidizer = app.UITable_R.Data{app.ind_Oxidizer, 2};
    else
        app.PD.S_Oxidizer = [];
    end
    % INERT
    if sum(app.ind_Inert)
        app.PD.S_Inert = regexprep(app.UITable_R.Data(app.ind_Inert,1),'\s+','')';
        app.PD.N_Inert = app.UITable_R.Data{app.ind_Inert, 2};
    else
        app.PD.S_Inert = [];
    end
end