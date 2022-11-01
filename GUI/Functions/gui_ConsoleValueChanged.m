function gui_ConsoleValueChanged(app, event)
    % Print output of commands through GUI's command window
    
    % Definitions
    FLAG_ERROR = false;
    % Read command
    try
        switch lower(app.Console.Value{1,1})
            case 'clear'
                ClearButtonPushed(app, event);
                output = ' ';
            case {'docs', 'documentation'}
                docs_CT;
                output = 'Redirecting to CT documentation using the default browser...';
            case {'license'}
                run('uilicense.mlapp');
                output = 'Running uilicense...';
            case {'update'}
                check_update(app.UIFigure);
                output = 'Checking for updates of Combustion Toolbox...';
            case {'validations'}
                run('uivalidations.mlapp');
                output = 'Running uivalidations...';
            case {'version'}
                [temp_tag, temp_date] = get_combustion_toolbox_version();
                output = sprintf('Version: %s\nDate: %s', temp_tag, temp_date);
            case {'save', 'export'}
                gui_save_results(app, app.export_results.format);
                output = 'Exporting results...';
            case {'settings', 'configuration'}
                uipreferences(app);
                output = 'Opening uipreferences...';
            case {'web', 'website', 'webpage'}
                website_CT;
                output = 'Redirecting to CT website using the default browser...';
            otherwise
                app.Console_text.Value = sprintf('Running: %s...', app.Console.Value{1,1});
                app.Lamp.Color = app.color_lamp_working;
                pause(0.1);
                output = evalc(app.Console.Value{1,1});
        end
    catch ME
        type = 'Error';
        output = sprintf('%s in function %s() at line %d.\nIdentifier: %s\nMessage: %s', ...
            type, ME.stack(1).name, ME.stack(1).line, ME.identifier, ME.message);
        app.Lamp.Color = app.color_lamp_error;
        FLAG_ERROR = true;
    end

    % Print output
    if ~isempty(output) || strcmp('clc', app.Console.Value{1,1})
        app.Console_text.Value = output;
    end
    app.Console.Value = '';

    % Set color lamp
    if ~FLAG_ERROR
        app.Lamp.Color = app.color_lamp_nothing;
    end
end