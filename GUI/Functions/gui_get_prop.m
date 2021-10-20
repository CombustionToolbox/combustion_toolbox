function [prop, flag_prop] = gui_get_prop(app, prop, name)
    if prop(1) == '['
        flag_prop = true;
        prop = sscanf(prop, '[%f:%f:%f]');
        prop = prop(1):prop(2):prop(3);
        app.PD.phi.value = app.PD.phi.value(1)*ones(1,numel(prop));
    else
        flag_prop = false;
        prop = sscanf(prop,'%f');
        if isempty(prop); prop = app.PD.(name).value; return; end
        app.PD.(name).value = prop;
    end
end