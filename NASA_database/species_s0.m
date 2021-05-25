function s0i = species_s0(Species,T,strThProp)
% global strThProp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates the entropy (in kJ/(mol·K)) at the specified temperature 
% T for the chemical species included in the following list:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% s0i = interp1(strThProp.(Species).T,strThProp.(Species).s0,T,'pchic');
% s0i = spline(strThProp.(Species).T,strThProp.(Species).s0,T);
% aux = griddedInterpolant(strThProp.(Species).T,strThProp.(Species).s0);
s0i = strThProp.(Species).s0curve(T)/1000;