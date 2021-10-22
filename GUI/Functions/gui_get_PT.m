function self = gui_get_PT(self)
    switch self.ProblemType.Value
        case '1' % * TP: Equilibrium composition at defined T and p
            self.PD.ProblemType = 'TP';
        case '2' % * HP: Adiabatic T and composition at constant p
            self.PD.ProblemType = 'HP';
        case '3' % * SP: Isentropic (i.e., adiabatic) compression/expansion to a specified p
            self.PD.ProblemType = 'SP';
        case '4' % * TV: Equilibrium composition at defined T and constant v
            self.PD.ProblemType = 'TV';  
        case '5' % * EV: Equilibrium composition at Adiabatic T and constant v
            self.PD.ProblemType = 'EV';
        case '6' % * SV: Isentropic (i.e., fast adiabatic) compression/expansion to a specified v
            self.PD.ProblemType = 'SV';
        case '7' % * SHOCK_I: CALCULATE PLANAR INCIDENT SHOCK WAVE
            self.PD.ProblemType = 'SHOCK_I';
        case '8' % * SHOCK_R: CALCULATE PLANAR POST-REFLECTED SHOCK STATE
            self.PD.ProblemType = 'SHOCK_R';
        case '9' % * DET: CALCULATE CHAPMAN-JOUGET STATE (CJ UPPER STATE)
            self.PD.ProblemType = 'DET';
        case '10' % * DET: CALCULATE OVERDRIVEN CHAPMAN-JOUGET STATE (CJ UPPER STATE)
            self.PD.ProblemType = 'DET_OVERDRIVEN';
    end
end