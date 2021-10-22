function self = SolveProblem(self, ProblemType)
% Save Problem Type
self.PD.ProblemType = ProblemType;
% Check inputs
self = check_inputs(self);
% Get Flags and length of the loop
self = get_FLAG_N(self);
self.C.l_phi = length(self.PD.phi.value);
for i=self.C.l_phi:-1:1
    % COMPUTE PROPERTIES INITIAL MIXTURE
    self = Define_FOI(self, i);
    % SOLVE SELECTED PROBLEM
    self = selectProblem(self, i);
    % DISPLAY RESULTS COMMAND WINDOW
    results(self, i);
end
end

% SUB-PASS FUNCTIONS
function self = selectProblem(self, i)
    % Solve selected problem
    switch self.PD.ProblemType
        case {'SHOCK_I', 'SHOCK_R'}
            u1 = self.PD.u1.value(i);
            if strcmp(self.PD.ProblemType,'SHOCK_I')
                if i==self.C.l_phi
                    [self.PS.strR{i}, self.PS.strP{i}] = shock_incident(self, self.PS.strR{i}, u1);
                else
                    [self.PS.strR{i}, self.PS.strP{i}] = shock_incident(self, self.PS.strR{i}, u1, self.PS.strP{i+1});
                end
            else
                if i==self.C.l_phi
                    [self.PS.strR{i}, self.PS.str2{i}, self.PS.strP{i}] = shock_reflected(self, self.PS.strR{i}, u1);
                else
                    [self.PS.strR{i}, self.PS.str2{i}, self.PS.strP{i}] = shock_reflected(self, self.PS.strR{i}, u1, self.PS.str2{i+1}, self.PS.strP{i+1});
                end
            end
        case {'DET'}
            if i==self.C.l_phi
                [self.PS.strR{i}, self.PS.strP{i}] = cj_detonation(self, self.PS.strR{i});
            else
                [self.PS.strR{i}, self.PS.strP{i}] = cj_detonation(self, self.PS.strR{i}, self.PS.strP{i+1});
            end
        case 'DET_OVERDRIVEN'
            overdriven = self.PD.overdriven.value(i);
            [self.PS.strR{i}, self.PS.strP{i}] = overdriven_detonation(self, self.PS.strR{i}, overdriven);
        otherwise
            if length(self.PD.pP.value) > 1
                pP = self.PD.pP.value(i);
            else
                pP = self.PD.pP.value;
            end
            if i == self.C.l_phi
                self.PS.strP{i} = equilibrate(self, self.PS.strR(i), pP);
            else
                self.PS.strP{i} = equilibrate(self, self.PS.strR(i), pP, self.PS.strP(i + 1));
            end
    end
end