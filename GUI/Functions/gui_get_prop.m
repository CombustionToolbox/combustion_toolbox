function [value, flag_prop] = gui_get_prop(app, name, value)
    if value(1) == '['
        flag_prop = true;
        value = sscanf(value, '[%f:%f:%f]');
        value = value(1):value(2):value(3);
        app.PD.phi.value = app.PD.phi.value(1)*ones(1,length(value));
    else
        flag_prop = false;
        value = sscanf(value,'%f');
        if isempty(value); value = app.PD.(name).value; return; end
        app.PD.(name).value = value;
    end
end