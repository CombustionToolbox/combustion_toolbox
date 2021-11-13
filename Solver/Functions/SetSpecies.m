function M = SetSpecies(self, species, N, T)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   function M = SetSpecies(self, species, N, T)
%   DESCRIPTION: Create the stoichiometric matrix
%   INPUT:
%       self    = contains necessary data
%       Species = species contained in the system 
%       N       = number of moles of each species
%       T       = Temperature of the species  
%   OUTPUT:
%       M = properties matrix           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = self.C.M0.value;
R0 = self.C.R0;

for n = length(N):-1:1
    mi = N(n) * self.DB.(species{n}).mm;
    mmi = self.DB.(species{n}).mm;
    hfi = self.DB.(species{n}).hf/1000;
    efi = self.DB.(species{n}).ef/1000;
    swtCondensed = self.DB.(species{n}).swtCondensed;
    if length(self.DB.(species{n}).T) > 1
        DhTi = species_DhT(species{n},T,self.DB);
        DeTi = species_DeT(species{n},T,self.DB);
        cPi  = species_cP(species{n},T,self.DB);
        cVi  = species_cV(species{n},T,self.DB);
        s0i  = species_s0(species{n},T,self.DB);
        if ~swtCondensed
            pVi = N(n)*R0*T/100; % For ideal gases
        else
            pVi = 0; % For condensed species
        end
    else
        DhTi = 0;
        DeTi = 0;
        cPi  = 0;
        cVi  = 0;
        s0i  = 0;
        pVi  = 0;
    end   
    ind = find_ind(species(n), self.S.LS);
    M(ind, :) = [N(n), N(n) * [hfi, DhTi, efi, DeTi, cPi, cVi, s0i], pVi, swtCondensed, mi, mmi];
end
