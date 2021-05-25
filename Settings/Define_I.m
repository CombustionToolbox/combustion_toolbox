function app = Define_I(app)
if ~isempty(app.PD.S_Inert)
    app.PD.R_Inert = SetSpecies(app.C.M0.value,app.PD.S_Inert,app.PD.N_Inert,app.PD.TR.value,find_ind(app.PD.S_Inert,app.S.LS), app.strThProp);
else
    app.PD.R_Inert = 0; % case without inert gases
end
if ~isempty(app.PD.S_Oxidizer) || ~isempty(app.PD.S_Inert)
    app.PS.strR_Oxidizer = ComputeProperties(app.C.A0.value,app.PD.R_Oxidizer+app.PD.R_Inert,app.PD.pR.value,app.PD.TR.value,app.E.ind_C,app.E.ind_H);
end