function gui_update_phi(app, self)
    % Update GUI: equivalence ratio, O/F, and percentage Fuel
    %
    % Args:
    %     app (object): Combustion Toolbox app object
    %     self (object): Data of the mixture, conditions, and databases
    
    try
        app.edit_phi.Value = sprintf('%.5g', round(self.PS.strR{1}.phi, 5));
    catch
        app.edit_phi.Value = self.PS.strR{1}.phi; 
    end

    app.edit_phi2.Value = app.edit_phi.Value;
    app.edit_OF.Value = self.PS.strR{1}.OF;
    app.edit_F.Value = self.PS.strR{1}.percentage_Fuel;
end