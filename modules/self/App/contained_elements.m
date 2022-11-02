function self = contained_elements(self)
    % Obtain containted elements from the given set of species (reactants and products)
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %
    % Returns:
    %     self (struct): Data of the mixture, conditions, and databases

    L_formula = [];

    for k = self.S.NS:-1:1
        L_E1 = []; L_E2 = [];
        formula = self.S.LS_formula{k};

        [idx0, idxf] = regexp(formula, "[A-Z]{2,}");

        for j = length(idxf):-1:1
            L_E2{j} = formula(idx0(j):idxf(j));
            formula(idx0(j):idxf(j)) = ' ';

            if isempty(L_E2)
                L_E2{j} = [];
            end

        end

        [~, idxf] = regexp(formula, "[A-Z]{1}");

        for j = length(idxf):-1:1
            L_E1{j} = formula(idxf(j));
        end

        L_formula = [L_formula, L_E1, L_E2];
    end

    self.E.elements = unique(L_formula);
    self.E.NE = numel(self.E.elements);
    self.E.ind_C = find(ismember(self.E.elements, 'C'));
    self.E.ind_H = find(ismember(self.E.elements, 'H'));
    self.E.ind_O = find(ismember(self.E.elements, 'O'));
    self.E.ind_N = find(ismember(self.E.elements, 'N'));
    self.E.ind_E = find(ismember(self.E.elements, 'E'));
    self.E.ind_S = find(ismember(self.E.elements, 'S'));
    self.E.ind_Si = find(ismember(self.E.elements, 'SI'));
end
