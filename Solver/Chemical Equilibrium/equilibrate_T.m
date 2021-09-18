function strP = equilibrate_T(self, strR, pP, TP)
    % Compute number of moles 
    [N, DeltaNP] = equilibrium(self, pP, TP, strR);
    % Compute thermodynamic derivates
    [dNi_T, dN_T] = equilibrium_dT(self, N, TP, strR);
    self.dNi_T = dNi_T;
    self.dN_T = dN_T;
    % Compute properties of all species
    P = SetSpecies(self, self.S.LS, N(:, 1), TP);

    if strfind(self.PD.ProblemType, 'P') == 2
        strP = ComputeProperties(self, P, pP, TP);
    else
        NP = sum(P(:, 1) .* (1 - P(:, 10)));
        pP = compute_pressure(self, strR, TP, NP);
        strP = ComputeProperties(self, P, pP, TP);
    end    
    strP.error_moles = DeltaNP;
end

function pP = compute_pressure(self, strR, TP, N)
    pP = (N * TP * self.C.R0 / (strR.v/1e3)) / 1e5;
end