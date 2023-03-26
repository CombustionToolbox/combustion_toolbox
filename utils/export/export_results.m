function export_results(self)
    % Export results of reactants (mix1) and products (mix2) into a .xls
    % file
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    
    % Check export_results value
    if ~self.Misc.export_results.value
        return
    end
    
    % Format data
    data_mix1 = FormattedOutput_test([], self.PD.phi.value, self.PS.strR, self.S.LS);
    data_mix2 = FormattedOutput_test([], self.PD.phi.value, self.PS.strP, self.S.LS);
    
    % Export data
    switch lower(self.Misc.export_results.format)
        case {'.xls', 'xls', 'excel', 'spreadsheet'}
            sheet_name = self.PD.ProblemType;
            filename = strcat(self.Misc.export_results.filename, '.xls');
            writecell(data_mix1, filename, 'Sheet', strcat('R-', sheet_name));
            writecell(data_mix2, filename, 'Sheet', strcat('P-', sheet_name));
        case {'.mat', 'mat', 'matlab'}
            filename = strcat(self.Misc.export_results.filename, '.mat');
            save(filename, 'data_mix1', 'data_mix2');
    end
    
end