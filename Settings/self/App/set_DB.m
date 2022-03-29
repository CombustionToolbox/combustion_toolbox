function self = set_DB(self, FLAG_REDUCED_DB, FLAG_FAST)
    % Generate Database with custom polynomials from DB_master
    %
    % Args:
    %     self (struct):          Data of the mixture, conditions, and databases    
    %     FLAG_REDUCED_DB (bool): Flag compute from reduced database
    %     FLAG_FAST (bool):       Flag load databases
    % 
    % Returns:
    %     self (struct):          Data of the mixture, conditions, and databases

    if ~FLAG_FAST
        % FLAG_REDUCED_DB:
        %   * False: complete Database
        %   * True:  reduced Database
        self.DB_master = generate_DB_master(FLAG_REDUCED_DB); 
        % Reduced Database with tabulated data of selected species
        self.DB = generate_DB(self.DB_master);
    end
    self.S.LS_DB = fieldnames(self.DB); 
    self.S.NS_DB = numel(self.S.LS_DB);
end