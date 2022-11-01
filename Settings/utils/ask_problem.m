function PT = ask_problem(self)
    % Create a list selection dialog box
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %
    % Returns:
    %     PT (string): String with the problem selected

    if self.Misc.FLAG_FOI && ~self.Misc.FLAG_GUI

        try
            fn = {'TP', 'HP', 'SP', 'TV', 'EV', 'SV', 'SHOCK_I', 'SHOCK_R', 'DET', 'DET_OVERDRIVEN'};
            [indx, ~] = listdlg('PromptString', 'Select a problem:', ...
                'SelectionMode', 'single', ...
                'ListString', fn, 'ListSize', [150, 150]);
            PT = fn{indx};
        catch
            error('Problem type not selected.')
        end

    else

        if self.Misc.FLAG_GUI
            self = gui_get_PT(self);
        end

        PT = self.PD.ProblemType;
    end

end
