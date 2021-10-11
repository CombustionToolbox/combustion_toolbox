function self = set_Database(self, reducedDB)
    % False: complete DataBase; True: reduced DB
    self.DB_master = ParseThermoInp(reducedDB); 
    % Reduced Database with tabulated data of selected species
    self.DB = GenerateDatabase(self.DB_master);
    self.S.LS_DB = fieldnames(self.DB); 
    self.S.NS_DB = numel(self.S.LS_DB);
end