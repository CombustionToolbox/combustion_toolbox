function g0i = species_g0_new(Species,T,strThProp)
% global strThProp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates the gibbs free energy (in kj/mol) at the specified temperature 
% T for the chemical species included in the following list:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% g0i = interp1(strThProp.(Species).T,strThProp.(Species).g0,T,'pchic');
% g0i = spline(strThProp.(Species).T,strThProp.(Species).g0,T);
% aux = griddedInterpolant(strThProp.(Species).T,strThProp.(Species).g0);
g0i = strThProp.(Species).g0curve(T)/1000;
%  tic, for T=300:100:3000 ; species_g0_new(strThProp,'CO',1800), end, toc