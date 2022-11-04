function gui_update_terminal(app, self, type)
    % Update GUI command window
    switch lower(type)
        case {'start'}
            message = sprintf('Solving %s problem... ', self.PD.ProblemType);
        case {'finish'}
            message = strcat(app.Console_text.Value, sprintf('Done! check tab "Results"\n\nElapsed time is %.6g seconds', self.Misc.timer_loop));
    end
    app.Console_text.Value = message;
end