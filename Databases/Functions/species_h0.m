function h0i = species_h0(species, T, DB)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates the enthalpy (in kJ/mol) at the specified temperature 
% T for the chemical species
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
    h0i = DB.(species).h0curve(T)/1000;
catch
    h0i = DB.(species).h0/1000;
end