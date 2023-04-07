function gui_update_terminal(app, self, type)
    % Update GUI command window depending on the type of message
    %
    % Args:
    %     app (object): Combustion Toolbox app object
    %     self (struct): Data of the mixture, conditions, and databases
    %     type (char): Type of message to be displayed

    switch lower(type)
        case {'start'}
            message = sprintf('Solving %s problem... ', self.PD.ProblemType);
        case {'finish'}
            message = strcat(app.Console_text.Value, sprintf('Done! check tab "Results"\n\nElapsed time is %.6g seconds', self.Misc.timer_loop));
    end

    app.Console_text.Value = message;
end