function cPi = species_cP_new(Species,T,strThProp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates the specific heat at constant pressure (in j/(mol·k)) for the 
% specified chemical species (Species) at the specified temperature (T)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cPi = interp1(strThProp.(Species).T,strThProp.(Species).cp,T,'pchic');
% cPi = spline(strThProp.(Species).T,strThProp.(Species).cp,T);
% aux = griddedInterpolant(strThProp.(Species).T,strThProp.(Species).cp);
% cPi = splinterp1(strThProp.(Species).cPcurve,T);
cPi = strThProp.(Species).cPcurve(T);
