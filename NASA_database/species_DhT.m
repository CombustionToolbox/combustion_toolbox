function DhTi = species_DhT(Species,T,strThProp)
% global strThProp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates the thermal enthalpy (in kj/mol) for the specified 
% chemical species (Species) at the specified temperature (T)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DhTi = interp1(strThProp.(Species).T,strThProp.(Species).DhT,T,'pchic')/1000;
% DhTi = spline(strThProp.(Species).T,strThProp.(Species).DhT,T)/1000;
% aux = griddedInterpolant(strThProp.(Species).T,strThProp.(Species).DhT);
DhTi = strThProp.(Species).DhTcurve(T)/1000;