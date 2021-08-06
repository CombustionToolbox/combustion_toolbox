function DhTi = species_DhT(Species,T,strThProp)
% global strThProp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates the thermal enthalpy (in kJ/mol) for the specified 
% chemical species (Species) at the specified temperature (T)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DhTi = strThProp.(Species).DhTcurve(T)/1000;