function gui_ConsoleValueChanged(app, event)
    if strcmp('clear', app.Console.Value{1,1})
        ClearButtonPushed(app, event);
        output = ' ';
    else
        try
            output = evalc(app.Console.Value{1,1});
        catch error
            output = error.message;
        end
    end
    if ~isempty(output) || strcmp('clc', app.Console.Value{1,1})
        app.Console_text.Value = output;
    end
    app.Console.Value = '';
end