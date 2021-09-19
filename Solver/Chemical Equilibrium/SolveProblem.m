function self = SolveProblem(self, i)

switch self.PD.ProblemType
    case {'SHOCK_I', 'SHOCK_R'}
        u1 = self.PD.u1.value(i);
        if strcmp(self.PD.ProblemType,'SHOCK_I')
            if i==self.C.l_phi
                [self.PS.strR{i}, self.PS.strP{i}] = SolveProblemSHOCK_I(self, self.PS.strR{i}, self.PD.pR.value, self.PD.TR.value, u1);
%                 [self.PS.strR{i}, self.PS.strP{i}] = shock_incident(self, self.PS.strR{i}, u1);
            else
                [self.PS.strR{i}, self.PS.strP{i}] = SolveProblemSHOCK_I_fast(self, self.PS.strR{i}, self.PD.pR.value, self.PD.TR.value, u1, self.PS.strP{i+1});
%                 [self.PS.strR{i}, self.PS.strP{i}] = shock_incident_fast(self, self.PS.strR{i}, u1, self.PS.strP{i+1});
            end
        else
            if i==self.C.l_phi
                [self.PS.strR{i}, self.PS.str2{i}, self.PS.strP{i}] = SolveProblemSHOCK_R(self, self.PS.strR{i}, self.PD.pR.value, self.PD.TR.value, u1);
            else
                [self.PS.strR{i}, self.PS.str2{i}, self.PS.strP{i}] = SolveProblemSHOCK_R_fast(self, self.PS.strR{i}, self.PD.pR.value, self.PD.TR.value, u1, self.PS.str2{i+1}, self.PS.strP{i+1});
            end
        end
    case {'DET'}
        [self.PS.strR{i},self.PS.strP{i},self.TN.guess] = SolveProblemDET(self, self.PS.strR{i},self.PD.pR.value,self.PD.TR.value, self.TN.guess);
    case 'DET_OVERDRIVEN'
        [self.PS.strR{i},self.PS.strP{i},self.TN.guess] = SolveProblemDET_OVERDRIVEN(self, self.PS.strR{i},self.PD.pR.value,self.PD.TR.value,self.TN.guess, self.PD.overdriven.value);
    otherwise
        if i == self.C.l_phi
            self.PS.strP{i} = equilibrate(self, self.PS.strR(i), self.PD.pP.value);
        else
            self.PS.strP{i} = equilibrate(self, self.PS.strR(i), self.PD.pP.value, self.PS.strP(i + 1));
        end
end
end