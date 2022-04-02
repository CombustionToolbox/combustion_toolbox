function s0i = species_s0(species, T, DB)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates the entropy (in kJ/(mol-K)) at the specified temperature 
% T for the chemical species included in the following list:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s0i = DB.(species).s0curve(T)/1000;