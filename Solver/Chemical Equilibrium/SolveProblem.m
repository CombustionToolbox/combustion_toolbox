function app = SolveProblem(app, i)
% Specify the problem type and the thermodynamic properties of the product
% mixture (P) required for the computations
switch app.PD.ProblemType
    case {'TP','TV'}
        app.PS.strP{i} = SolveProblemTP_TV(app, app.PS.strR{i},app.PD.phi.value(i),app.PD.pR.value, app.PD.TP.value);
    case {'HP'}
        if i==app.C.l_phi
            app.PS.strP{i} = SolveProblemHP(app, app.PS.strR{i}, app.PD.phi.value(i), app.PD.pR.value);
        else
            app.PS.strP{i} = SolveProblemHP_EV_fast(app, app.PS.strR{i},app.PD.phi.value(i),app.PD.pR.value,app.PS.strP{i+1});
        end
    case {'EV'}
        if i==app.C.l_phi
            app.PS.strP{i} = SolveProblemEV(app, app.PS.strR{i}, app.PD.phi.value(i), app.PD.pR.value);
        else
            app.PS.strP{i} = SolveProblemHP_EV_fast(app, app.PS.strR{i},app.PD.phi.value(i),app.PD.pR.value,app.PS.strP{i+1});
        end
    case 'SP'
        pP = app.PD.pP.value(i);
        if i==app.C.l_phi
            app.PS.strP{i} = SolveProblemSP(app, app.PS.strR{i}, app.PD.phi.value(i), pP);
        else
            app.PS.strP{i} = SolveProblemSP_fast(app, app.PS.strR{i}, app.PD.phi.value(i), pP, app.PS.strP{i+1});
        end
    case 'SV'
        vP_vR = app.PD.vP_vR.value(i);
        app.PS.strP{i} = SolveProblemSV(app, app.PS.strR{i}, app.PD.phi.value(i), vP_vR*app.PS.strR{i}.v);
    case {'SHOCK_I','SHOCK_R'}
        u1 = app.PD.u1.value(i);
        if strcmp(app.PD.ProblemType,'SHOCK_I')
            if i==app.C.l_phi
                [app.PS.strR{i}, app.PS.strP{i}] = SolveProblemSHOCK_I(app, app.PS.strR{i}, app.PD.phi.value(i), app.PD.pR.value, app.PD.TR.value, u1);
            else
                [app.PS.strR{i}, app.PS.strP{i}] = SolveProblemSHOCK_I_fast(app, app.PS.strR{i}, app.PD.phi.value(i), app.PD.pR.value, app.PD.TR.value, u1, app.PS.strP{i+1});
            end
        else
            if i==app.C.l_phi
                [app.PS.strR{i}, app.PS.str2{i}, app.PS.strP{i}] = SolveProblemSHOCK_R(app, app.PS.strR{i}, app.PD.phi.value(i), app.PD.pR.value, app.PD.TR.value, u1);
            else
                [app.PS.strR{i}, app.PS.str2{i}, app.PS.strP{i}] = SolveProblemSHOCK_R_fast(app, app.PS.strR{i}, app.PD.phi.value(i), app.PD.pR.value, app.PD.TR.value, u1, app.PS.str2{i+1}, app.PS.strP{i+1});
            end
        end
    case 'DET'
%         app.PD.TR.value = app.PD.TR_vector.value(i);
%         [app.PS.strR{i},app.PS.strP{i},app.TN.guess] = SolveProblemDET_main(app.PS.strR{i},app.PD.phi.value(i),app.PD.pR.value,app.PD.TR.value,app.TN.guess,app.E,app.S,app.C,app.M,app.PD,app.TN,app.strThProp);
        [app.PS.strR{i},app.PS.strP{i},app.TN.guess] = SolveProblemDET_main_2(app, app.PS.strR{i},app.PD.phi.value(i),app.PD.pR.value,app.PD.TR.value, app.TN.guess);
%         [app.PS.strR{i},app.PS.strP{i},app.TN.guess] = SolveProblemDET_hybrid_main(app.PS.strR{i},app.PD.phi.value(i),app.PD.pR.value,app.PD.TR.value,app.TN.guess,app.E,app.S,app.C,app.M,app.PD,app.TN,app.strThProp);
%         [app.PS.strR{i},app.PS.strP{i},app.TN.guess] = SolveProblemDET_Jouget_main(app.PS.strR{i},app.PD.phi.value(i),app.PD.pR.value,app.PD.TR.value,app.TN.guess,app.E,app.S,app.C,app.M,app.PD,app.TN,app.strThProp);
    case 'DET_OVERDRIVEN'
        [app.PS.strR{i},app.PS.strP{i},app.TN.guess] = SolveProblemDET_OVERDRIVEN(app, app.PS.strR{i},app.PD.phi.value(i),app.PD.pR.value,app.PD.TR.value,app.TN.guess, app.PD.overdriven.value);
end
end