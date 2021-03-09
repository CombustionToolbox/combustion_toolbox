function M = SetSpecies(M,S,N,T,id,strThProp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   function M = SetSpecies_new(A0,N,T)
%   DESCRIPTION: Create the stoichiometric matrix
%   INPUT:
%       A0 = unity stoichiometric matrix
%       S  = species contained in the system 
%       N  = number of moles of each specie  
%       T  = Temperature of the specie       
%       id = index of the species            
%   OUTPUT:
%       M = properties matrix           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R0 = 8.3144598; % [J/(K mol)]. Universal gas constant
M(id,1) = N;

for n = length(N):-1:1
    hfi = strThProp.(S{n}).hf/1000;
    efi = strThProp.(S{n}).ef/1000;
    if length(strThProp.(S{n}).T) > 1
        DhTi = species_DhT_new(S{n},T,strThProp);
        DeTi = species_DeT_new(S{n},T,strThProp);
        cPi = species_cP_new(S{n},T,strThProp);
        cVi = species_cV_new(S{n},T,strThProp);
        s0i = species_s0_new(S{n},T,strThProp);
        swtCondensed = strThProp.(S{n}).swtCondensed;
        mi = N(n)*strThProp.(S{n}).mm;
        mmi = strThProp.(S{n}).mm;
        if swtCondensed == 0
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
        swtCondensed = strThProp.(S{n}).swtCondensed;
        pVi = 0;
    end
    M(id(n),(2:end)) = [N(n)*[hfi, DhTi, efi, DeTi, cPi, cVi, s0i], pVi, swtCondensed, mi, mmi];
end
