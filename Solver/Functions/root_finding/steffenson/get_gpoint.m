function [gpoint, gpoint_relative]= get_gpoint(self, strR, pP, attr_name, x)
    strP = equilibrate_T(self, strR, pP, x);
    gpoint = strP.(attr_name) - strR.(attr_name);
    gpoint_relative = gpoint / strP.(attr_name);
end