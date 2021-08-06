function M = SetSpecies(self, Species, N, T)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   function M = SetSpecies_new(A0,N,T)
%   DESCRIPTION: Create the stoichiometric matrix
%   INPUT:
%       A0 = unity stoichiometric matrix
%       Species  = species contained in the system 
%       N  = number of moles of each specie  
%       T  = Temperature of the specie       
%       id = index of the species            
%   OUTPUT:
%       M = properties matrix           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = self.C.M0.value;
R0 = self.C.R0;

for n = length(N):-1:1
    hfi = self.strThProp.(Species{n}).hf/1000;
    efi = self.strThProp.(Species{n}).ef/1000;
    if length(self.strThProp.(Species{n}).T) > 1
        DhTi = species_DhT(Species{n},T,self.strThProp);
        DeTi = species_DeT(Species{n},T,self.strThProp);
        cPi = species_cP(Species{n},T,self.strThProp);
        cVi = species_cV(Species{n},T,self.strThProp);
        s0i = species_s0(Species{n},T,self.strThProp);
        swtCondensed = self.strThProp.(Species{n}).swtCondensed;
        mi = N(n)*self.strThProp.(Species{n}).mm;
        mmi = self.strThProp.(Species{n}).mm;
        if ~swtCondensed
            pVi = N(n)*R0*T/100; % For ideal gases
        else
            pVi = 0; % For condensed species
        end
    else
        DhTi = 0;
        DeTi = 0;
        cPi = 0;
        cVi = 0;
        s0i = 0;
        swtCondensed = self.strThProp.(Species{n}).swtCondensed;
        pVi = 0;
    end
    ind = find_ind(Species(n), self.S.LS);
    M(ind,:) = [N(n), N(n)*[hfi, DhTi, efi, DeTi, cPi, cVi, s0i], pVi, swtCondensed, mi, mmi];
end
