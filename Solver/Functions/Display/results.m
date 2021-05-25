function results(app, i)
    if ~strcmp(app.PD.ProblemType,'SHOCK_R') && ~strcmp(app.PD.ProblemType,'DET_OVERDRIVEN')
        displayresults(app.PS.strR{i},app.PS.strP{i},app.PD.ProblemType,app.C.mintol_display,app.S.LS);
    elseif ~strcmp(app.PD.ProblemType,'DET_OVERDRIVEN')
        displayresults(app.PS.strR{i},app.PS.str2{i},app.PS.strP{i},app.PD.ProblemType,app.C.mintol_display,app.S.LS); % Display all results SHOCK_R
    end
end