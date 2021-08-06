function self = set_Database(self, reducedDB)
    % False: complete DataBase; True: reduced DB
    self.strMaster = ParseThermoInp(reducedDB); 
    % Struct with tabulated data of selected species
    self.strThProp = GenerateDatabase(self.strMaster);
    self.S.LS_DB = fieldnames(self.strThProp); 
    self.S.NS_DB = numel(self.S.LS_DB);
end