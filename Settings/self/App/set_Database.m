function self = set_Database(self, reducedDB)
    % False: complete DataBase; True: reduced DB
    self.strMaster = ParseThermoInp(reducedDB); 
    % Struct with tabulated data of selected species
    self.strThProp = GenerateDatabase(self.strMaster);
    self.S.LS_DB = fieldnames(self.strThProp); 
    self.S.NS_DB = numel(self.S.LS_DB);
    self.S.formula_DB = get_formula(self.S.LS_DB, self.strThProp);
end

% NESTED FUNCTIONS
function formula_DB = get_formula(LS_DB, DB)
    for i=length(LS_DB):-1:1
        formula_DB{i} = DB.(LS_DB{i}).txFormula;
    end
end