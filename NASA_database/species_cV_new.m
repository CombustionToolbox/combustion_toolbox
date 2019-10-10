function cVi = species_cV_new(Species,T,strThProp)
% global strThProp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates the specific heat at constant volume (in j/(mol·k)) for the 
% specified chemical species (Species) at the specified temperature (T)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cVi = interp1(strThProp.(Species).T,strThProp.(Species).cv,T,'pchic');
% cVi = spline(strThProp.(Species).T,strThProp.(Species).cv,T);
% aux = griddedInterpolant(strThProp.(Species).T,strThProp.(Species).cv);
cVi = strThProp.(Species).cVcurve(T);