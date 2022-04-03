% In equilibrium_dT

dpi_T = x(temp_NS+1:end-1);
for i=temp_NE:-1:1
   cP_prime_1(i) = sum(A0(temp_ind_nswt, temp_ind_E(i)) .* moles(temp_ind_nswt, 1) .* h0(temp_ind_nswt));
end

cP_prime = (sum(cP_prime_1 .* dpi_T') + sum(h0(temp_ind_swt) .* dNi_T(temp_ind_swt)) +...
    sum(h0(temp_ind_nswt) .* moles(temp_ind_nswt, 1)) * dN_T + sum(moles(temp_ind_nswt, 1) .* h0(temp_ind_nswt).^2 / (C.R0 * T))) / T;

% In ComputeProperties
mix.cP = self.dpi_T + sum(SpeciesMatrix(:,6));