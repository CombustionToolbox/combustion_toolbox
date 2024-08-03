function gui_update_terminal(app, type, varargin)
    % Update GUI command window depending on the type of message
    %
    % Args:
    %     app (object): Combustion Toolbox app object
    %     type (char): Type of message to be displayed
    %
    % Additional inputs:
    %     * elapsedTime (float): Time expended in the solver

    switch lower(type)
        case {'start'}
            message = sprintf('Solving %s problem... ',app.ProblemType.Value);
        case {'finish'}
            elapsedTime = varargin{1};
            message = sprintf('%sDone! check tab "Results"\n\nElapsed time is %.6g seconds', app.Console_text.Value{:}, elapsedTime);
    end
    
    % Update GUI command window message
    app.Console_text.Value = message;
end