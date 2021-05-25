function app = Define_O(app)
if ~isempty(app.PD.S_Oxidizer)
    app.PD.R_Oxidizer = SetSpecies(app.C.M0.value,app.PD.S_Oxidizer,app.PD.N_Oxidizer,app.PD.TR.value,find_ind(app.PD.S_Oxidizer,app.S.LS), app.strThProp);
else
    app.PD.R_Oxidizer = 0;
end
