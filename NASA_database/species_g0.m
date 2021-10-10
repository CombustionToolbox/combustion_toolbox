function g0i = species_g0(species, T, DB)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates the gibbs free energy (in kJ/mol) at the specified temperature 
% T for the chemical species
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
    g0i = DB.(species).g0curve(T)/1000;
catch
    g0i = DB.(species).g0/1000; % cryogenic species
end