function self = Check_FOI(self, FOI_species)
    if self.Misc.FLAG_FOI
        self.Misc.FLAG_FOI = false;
        FLAG = false;
        for i=1:numel(FOI_species)
            if ~strcmp(self.S.LS, FOI_species(i))                 
                self.S.LS = [self.S.LS, FOI_species(i)];
                try 
                    self.M.ind_minor = [self.M.ind_minor, self.M.ind_minor(end) + 1];
                catch
                    self.M.ind_minor = 1;
                end
                self.M.Lminors = numel(self.S.LS);
                FLAG = true;
            end
        end
        if FLAG
            self.S.ind_nswt = [];
            self.S.ind_swt = [];
            self = Initialize(self);
        end
        self.Misc.FLAG_FIRST = false;
    end
end