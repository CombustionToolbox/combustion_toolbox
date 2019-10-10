function [t,C,Xi] = NO_Zeldovich(T,p,t0,tf,C0)
R = 8.314472;            % Ideal gas constant [J/mol K]
A = [1.8e14,6.4e9,3e13]; % Frecuency factor [1/s]
Ea= 1e3*[-319,-26.1];    % Activation energy [kcal/mol]
[t,C]=ode23t(@(t,y)NO_mechanism(t,y,T,A,Ea,R),[t0 tf],C0);
for i=length(C):-1:1
    C_total(i) = sum(C(i,:));
end
Xi = C./C_total';
t = t./t(end);