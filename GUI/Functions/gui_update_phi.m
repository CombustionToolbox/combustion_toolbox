function gui_update_phi(obj, app)
    % Update GUI: equivalence ratio, O/F, and percentage Fuel
    try
        obj.edit_phi.Value = sprintf('%.5g', round(app.PS.strR{1}.phi, 5));
    catch
        obj.edit_phi.Value = app.PS.strR{1}.phi; 
    end
    obj.edit_phi2.Value = obj.edit_phi.Value;
    obj.edit_OF.Value = 1/app.PS.strR{1}.FO;
    obj.edit_F.Value = app.PS.strR{1}.percentage_Fuel;
end