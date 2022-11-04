function app = gui_empty_UITables(app)
    % Clear data UITables and set to default the value of the equivalence ratio (-)
    app.UITable_R.Data  = [];
    app.UITable_P.Data  = [];
    app.UITable_R2.Data = [];
    app.edit_phi.Value = '-';
end