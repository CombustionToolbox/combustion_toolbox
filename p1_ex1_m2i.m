app.PD.N_Fuel = 1; % moles of Fuel
app.M.minor_products = {'O','OH','H'};
app.C.l_phi = length(app.PD.phi.Value);

for i=app.C.l_phi:-1:1
    app.PD.R_Fuel = SetSpecies(app.C.M0.Value,app.PD.S_Fuel,app.PD.N_Fuel,app.PD.TR.Value,find_idx(app.PD.S_Fuel,app.S.NameSpecies),strThProp);
    app.PS.strR_Fuel{i} = ComputeProperties(app.C.A0.Value,app.PD.R_Fuel,app.PD.pR.Value,app.PD.TR.Value,app.E.ind_C,app.E.ind_H); app.PS.strR{i}.phi = app.PD.phi.Value(i);
    x = app.PS.strR_Fuel{i}.NatomE(app.E.ind_C);
    y = app.PS.strR_Fuel{i}.NatomE(app.E.ind_H);
    z = app.PS.strR_Fuel{i}.NatomE(app.E.ind_O);
    app.PD.phi_t = x+y/4-z/2;
    app.PD.S_Oxidizer = {'O2'}; app.PD.N_Oxidizer = app.PD.phi_t/app.PD.phi.Value(i);
    app.PD.R_Oxidizer = SetSpecies(app.C.M0.Value,app.PD.S_Oxidizer,app.PD.N_Oxidizer,app.PD.TR.Value,find_idx(app.PD.S_Oxidizer,app.S.NameSpecies),strThProp);
    app.PD.S_Inert = {'N2'}; app.PD.N_Inert = app.PD.phi_t/app.PD.phi.Value(i)*79/21;
    app.PD.R_Inert = SetSpecies(app.C.M0.Value,app.PD.S_Inert,app.PD.N_Inert,app.PD.TR.Value,find_idx(app.PD.S_Inert,app.S.NameSpecies),strThProp);
    R = app.PD.R_Fuel + app.PD.R_Oxidizer + app.PD.R_Inert;
    app.PS.strR{i} = ComputeProperties(app.C.A0.Value,R,app.PD.pR.Value,app.PD.TR.Value,app.E.ind_C,app.E.ind_H); app.PS.strR{i}.phi = app.PD.phi.Value(i);
    if i==app.C.l_phi
        app.PS.strP{i} = SolveProblemHP_EV(app.PS.strR{i},app.PD.phi.Value(i),app.PD.pR.Value,app.E,app.S,app.C,app.M,app.PD,app.TN,strThProp);
    else
        app.PS.strP{i} = SolveProblemHP_EV_fast(app.PS.strR{i},app.PD.phi.Value(i),app.PD.pR.Value,app.PS.strP{i+1},app.E,app.S,app.C,app.M,app.PD,app.TN,strThProp);
    end
    T(i) = app.PS.strP{i}.T;
    displayresults_m2i(app.PS.strR{i},app.PS.strP{i},app.PD.ProblemType,app.C.mintol_display,app.S.NameSpecies);
end
fprintf('Stoichiometric equivalence ratio = %.2f\n',app.PD.phi_t);
displaysweepresults(app.PS.strR,app.PD.phi.Value,app.S.NameSpecies,app.C.mintol_display);
