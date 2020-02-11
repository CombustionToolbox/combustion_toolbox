function W = Compute_W(strR,strR_Fuel)
for i=length(strR):-1:1
    Y_Fuel = strR_Fuel.mi/strR{i}.mi;
    W_Air  = strR{i}.mi-strR_Fuel.mi;
    W_Fuel = strR_Fuel.mi;
    W(i)  = (1-(W_Air/W_Fuel))/(1-(1-W_Air/W_Fuel)*Y_Fuel);    
end
