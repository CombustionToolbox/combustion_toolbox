function g0i = species_g0(Species,T,strThProp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates the gibbs free energy (in J/mol) at the specified temperature 
% T for the chemical species included in the following list:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g0i = strThProp.(Species).g0curve(T);