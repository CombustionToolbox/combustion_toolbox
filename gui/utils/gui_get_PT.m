function app = gui_get_PT(app)
    switch app.ProblemType.Value
        case '1' % * TP: Equilibrium composition at defined T and p
            app.PD.ProblemType = 'TP';
        case '2' % * HP: Adiabatic T and composition at constant p
            app.PD.ProblemType = 'HP';
        case '3' % * SP: Isentropic (i.e., adiabatic) compression/expansion to a specified p
            app.PD.ProblemType = 'SP';
        case '4' % * TV: Equilibrium composition at defined T and constant v
            app.PD.ProblemType = 'TV';  
        case '5' % * EV: Equilibrium composition at Adiabatic T and constant v
            app.PD.ProblemType = 'EV';
        case '6' % * SV: Isentropic (i.e., fast adiabatic) compression/expansion to a specified v
            app.PD.ProblemType = 'SV';
        case '7' % * SHOCK_I: CALCULATE PLANAR INCIDENT SHOCK WAVE
            app.PD.ProblemType = 'SHOCK_I';
        case '8' % * SHOCK_R: CALCULATE PLANAR POST-REFLECTED SHOCK STATE
            app.PD.ProblemType = 'SHOCK_R';
        case '9' % * DET: CALCULATE CHAPMAN-JOUGET STATE (CJ UPPER STATE)
            app.PD.ProblemType = 'DET';
        case '10' % * DET: CALCULATE OVERDRIVEN CHAPMAN-JOUGET STATE (CJ UPPER STATE)
            app.PD.ProblemType = 'DET_OVERDRIVEN';
    end
end