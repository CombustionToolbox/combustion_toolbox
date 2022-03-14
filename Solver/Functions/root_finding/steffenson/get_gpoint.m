function [gpoint, gpoint_relative] = get_gpoint(self, strR, pP, attr_name, x)
    strP = equilibrate_T(self, strR, pP, x);
    gpoint = (strP.(attr_name) - strR.(attr_name));
    gpoint_relative = gpoint / (strP.(attr_name));
    if strcmpi(attr_name, 's')
        gpoint = gpoint * 1e3;
        gpoint_relative = gpoint_relative * 1e-3;
    end
end