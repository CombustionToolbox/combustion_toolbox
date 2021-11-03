function self = set_DB(self, reducedDB, FLAG_FAST)
    if ~FLAG_FAST
        % False: complete DataBase; True: reduced DB
        self.DB_master = generate_DB_master(reducedDB); 
        % Reduced Database with tabulated data of selected species
        self.DB = generate_DB(self.DB_master);
    end
    self.S.LS_DB = fieldnames(self.DB); 
    self.S.NS_DB = numel(self.S.LS_DB);
end