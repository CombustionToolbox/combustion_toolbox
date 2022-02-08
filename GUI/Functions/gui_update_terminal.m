function gui_update_terminal(obj, app, type)
    % Update GUI command window
    switch lower(type)
        case {'start'}
            message = sprintf('Solving %s problem... ', app.PD.ProblemType);
        case {'finish'}
            app.Misc.timer_loop = toc(app.Misc.timer_0);
            message = strcat(obj.Console_text.Value, sprintf('Done! check tab "Results"\n\nElapsed time is %.6g seconds', app.Misc.timer_loop));
    end
    obj.Console_text.Value = message;
end