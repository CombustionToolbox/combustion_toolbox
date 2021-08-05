function cPi = species_cP(Species,T,strThProp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates the specific heat at constant pressure (in J/(mol-K)) for the 
% specified chemical species (Species) at the specified temperature (T)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cPi = strThProp.(Species).cPcurve(T);
