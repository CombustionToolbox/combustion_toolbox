function results(self, i)
    % Display results in the command window
    if self.Misc.FLAG_RESULTS
        if ~strcmpi(self.PD.ProblemType, 'SHOCK_R') && ~strcmpi(self.PD.ProblemType, 'ROCKET') && ~strcmpi(self.PD.ProblemType, 'DET_R') 
            displayresults(self.PS.strR{i}, self.PS.strP{i}, self.PD.ProblemType, self.C.mintol_display, self.S.LS);
        elseif strcmpi(self.PD.ProblemType, 'ROCKET') 
            displayresults(self.PS.strR{i}, self.PS.str2{i}, self.PS.strP{i}, self.PD.ProblemType, self.C.mintol_display, self.S.LS);
        else
            % Display results reflected shocks (SHOCK_R or DET_R)
            displayresults(self.PS.strR{i}, self.PS.str2{i}, self.PS.strP{i}, self.PD.ProblemType, self.C.mintol_display, self.S.LS);
        end
    end
end