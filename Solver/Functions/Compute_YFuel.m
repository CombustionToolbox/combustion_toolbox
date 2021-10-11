function YFuel = Compute_YFuel(strR,strR_Fuel)
for i=length(strR):-1:1
    YFuel(i) = strR_Fuel.mi/strR{i}.mi;
end
