function gui_fontcolor_error(app, item, MAX_ERROR)
    % Change fontcolor of the given item if the value is greater than
    % error_max
    value = app.(item).Value;
    if value > MAX_ERROR
        app.(item).FontColor = [0.6350 0.0780 0.1840];
    else
        app.(item).FontColor = [0 0 0];
    end
end