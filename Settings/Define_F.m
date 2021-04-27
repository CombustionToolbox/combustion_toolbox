function app = Define_F(app)
if isfield(app.PD,'S_Fuel')
    app.PD.R_Fuel = SetSpecies(app.C.M0.Value,app.PD.S_Fuel,app.PD.N_Fuel,app.PD.TR.Value,find_idx(app.PD.S_Fuel,app.S.NameSpecies), app.strThProp);
    app.PS.strR_Fuel = ComputeProperties(app.C.A0.Value,app.PD.R_Fuel,app.PD.pR.Value,app.PD.TR.Value,app.E.ind_C,app.E.ind_H);
    app.PD.Fuel.x = app.PS.strR_Fuel.NatomE(app.E.ind_C); 
    app.PD.Fuel.y = app.PS.strR_Fuel.NatomE(app.E.ind_H); 
    app.PD.Fuel.z = app.PS.strR_Fuel.NatomE(app.E.ind_O); app.PS.strR_Fuel.z = app.PD.Fuel.z;
    app.PD.Fuel.w = app.PS.strR_Fuel.NatomE(app.E.ind_N); app.PS.strR_Fuel.w = app.PD.Fuel.w;
    app.PD.Fuel.eps = 0;
    app.PD.phi_t = app.PD.Fuel.x+app.PD.Fuel.y/4-app.PD.Fuel.z/2;
end

