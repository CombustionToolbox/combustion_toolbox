function self = gui_get_parameters(self)
    % GET EQUIVALENCE RATIO
    [self.PD.phi.value, self.flag_phi] = gui_get_prop(self, self.edit_phi.Value, 'phi');
    if isempty(self.PD.phi.value)
        self.PD.phi.value = 1;
    end
    % GET CONDITIONS
    [self.PR1_vector, self.flag_PR1] = gui_get_prop(self, self.PR1.Value, self.PR1_var_name);
    [self.PR2_vector, self.flag_PR2] = gui_get_prop(self, self.PR2.Value, self.PR2_var_name);
    [self.PP1_vector, self.flag_PP1] = gui_get_prop(self, self.PP1.Value, self.PP1_var_name);
    [self.PP2_vector, self.flag_PP2] = gui_get_prop(self, self.PP2.Value, self.PP2_var_name);
    if strcmpi(self.PD.ProblemType, 'DET_OVERDRIVEN')
        [self.PR3_vector, self.flag_PR3] = gui_get_prop(self, self.PR3.Value, self.PR3_var_name);
    end
end