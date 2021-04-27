function app = Define_O(app)
if isfield(app.PD,'S_Oxidizer')
    app.PD.R_Oxidizer = SetSpecies(app.C.M0.Value,app.PD.S_Oxidizer,app.PD.N_Oxidizer,app.PD.TR.Value,find_idx(app.PD.S_Oxidizer,app.S.NameSpecies), app.strThProp);
else
    app.PD.R_Oxidizer = 0;
end
