function gui_fontcolor_error(app, item, MAX_ERROR)
    % Change fontcolor of the given item if the value is greater than
    % error_max
    value = app.(item).Value;

    if value > MAX_ERROR
        app.(item).FontColor = app.color_lamp_error;
        app.Lamp.Color = app.color_lamp_error;
    else
        app.(item).FontColor = [0 0 0];
    end

end