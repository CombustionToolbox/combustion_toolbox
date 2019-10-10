function DeTi = species_DeT_new(Species,T,strThProp)
% global strThProp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates the thermal internal energy (in kj/mol) for the specified 
% chemical species (Species) at the specified temperature (T)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DeTi = interp1(strThProp.(Species).T,strThProp.(Species).DeT,T,'pchic')/1000;
% DeTi = spline(strThProp.(Species).T,strThProp.(Species).DeT,T)/1000;
% aux = griddedInterpolant(strThProp.(Species).T,strThProp.(Species).DeT);
DeTi = strThProp.(Species).DeTcurve(T)/1000;
