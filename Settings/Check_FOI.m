function app = Check_FOI(app, FOI_species)
    if app.Misc.FLAG_FOI
        app.Misc.FLAG_FOI = false;
        FLAG = false;
        for i=1:numel(FOI_species)
            if ~strcmp(app.S.LS, FOI_species(i))                 
                app.M.minors_products = [app.M.minors_products, FOI_species(i)];
                try 
                    app.M.ind_minor = [app.M.ind_minor, app.M.ind_minor(end) + 1];
                catch
                    app.M.ind_minor = 1;
                end
                app.M.Lminors = numel(app.M.minors_products);
                FLAG = true;
            end
        end
        if FLAG
            app.S.ind_nswt = [];
            app.S.ind_swt = [];
            app = Initialize(app);
        end
        app.Misc.FLAG_FIRST = false;
    end
end