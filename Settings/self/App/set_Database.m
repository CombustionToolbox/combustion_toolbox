function self = set_Database(self, reducedDB)
    % False: complete DataBase; True: reduced DB
    self.strMaster = ParseThermoInp(reducedDB); 
    % Struct with tabulated data of selected species
    self.strThProp = GenerateDatabase(self.strMaster);
    self.S.namespecies = fieldnames(self.strThProp); 
    self.S.Nspecies = numel(self.S.namespecies);
end