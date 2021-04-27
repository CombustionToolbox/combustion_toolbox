function app = Define_FOI(app, i)
% COMPUTE PROPERTIES OF THE REACTIVES FOR THE GIVEN CONDITIONS
R = app.PD.R_Fuel+app.PD.R_Oxidizer+app.PD.R_Inert;
app.PS.strR{i} = ComputeProperties(app.C.A0.Value,R,app.PD.pR.Value,app.PD.TR.Value,app.E.ind_C,app.E.ind_H); app.PS.strR{i}.phi = app.PD.phi.Value(i);