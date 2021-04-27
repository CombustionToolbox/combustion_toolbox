function app = Define_I(app)
if isfield(app.PD,'S_Inert')
    app.PD.R_Inert = SetSpecies(app.C.M0.Value,app.PD.S_Inert,app.PD.N_Inert,app.PD.TR.Value,find_idx(app.PD.S_Inert,app.S.NameSpecies), app.strThProp);
else
    app.PD.R_Inert = 0; % case without inert gases
end
if isfield(app.PD,'S_Inert') || isfield(app.PD,'S_Inert')
    app.PS.strR_Oxidizer = ComputeProperties(app.C.A0.Value,app.PD.R_Oxidizer+app.PD.R_Inert,app.PD.pR.Value,app.PD.TR.Value,app.E.ind_C,app.E.ind_H);
end