function obj = gui_empty_UITables(obj)
    % Clear data UITables and set to default the value of the equivalence ratio (-)
    obj.UITable_R.Data  = [];
    obj.UITable_P.Data  = [];
    obj.UITable_R2.Data = [];
    obj.edit_phi.Value = '-';
end