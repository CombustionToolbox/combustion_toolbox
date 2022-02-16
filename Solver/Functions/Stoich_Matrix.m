function self = Stoich_Matrix(self)
    % Create stoichiometric matrix
    self.C.A0.value = zeros(self.S.NS, self.E.NE);
    self.C.M0.value = zeros(self.S.NS, self.C.N_prop.value);
    for i=1:self.S.NS
        txFormula = self.DB.(self.S.LS{i}).txFormula;
        self.DB.(self.S.LS{i}).Element_matrix = set_element_matrix(txFormula,self.E.elements);
        self.C.A0.value(i,self.DB.(self.S.LS{i}).Element_matrix(1,:)) = self.DB.(self.S.LS{i}).Element_matrix(2,:);
        self.C.M0.value(i,10) = self.DB.(self.S.LS{i}).swtCondensed;    
    end
    self.C.N0.value = self.C.M0.value(:, [1, 10]);
end