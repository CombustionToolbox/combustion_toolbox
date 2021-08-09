function gpoint = get_gpoint(self, strR, pP, attr_name, x)
    strP = equilibrate_T(self, strR, pP, x);
    gpoint = getattr(strP, attr_name) - getattr(strR, attr_name);
end