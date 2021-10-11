function h0i = species_h0(Species,T,strThProp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates the enthalpy (in kJ/mol) at the specified temperature 
% T for the chemical species
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
    h0i = strThProp.(Species).h0curve(T)/1000;
catch
    h0i = strThProp.(Species).h0/1000;
end