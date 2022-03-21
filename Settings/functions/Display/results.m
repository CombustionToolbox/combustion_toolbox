function results(self, i)
    % Display results in the command window
    if self.Misc.FLAG_RESULTS
        if ~isfield(self.PS, 'str2')
            displayresults(self.PS.strR{i}, self.PS.strP{i}, self.PD.ProblemType, self.C.mintol_display, self.S.LS);
        else
            % Display results reflected shocks (SHOCK_R, DET_R, or DET_OVERDRIVEN_R)
            displayresults(self.PS.strR{i}, self.PS.str2{i}, self.PS.strP{i}, self.PD.ProblemType, self.C.mintol_display, self.S.LS);
        end
    end
end