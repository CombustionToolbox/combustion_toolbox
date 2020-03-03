% load DET_H2_fixed2;

[rho, ~, W] = plot_density_YFuel(app.PS.strR,app.PS.strR_Fuel);
[q, YFuel, H, YFuel_st] = plot_heatrelease_YFuel(app.PS.strR,app.PS.strP,app.PS.strR_Fuel);


fileID = fopen(strcat('C',sprintf('%d',app.PS.strR{1}.x),'H',sprintf('%d',app.PS.strR{1}.y),'.txt'),'w');
for i=1:length(rho)
    fprintf(fileID,'%f   %f   %f\n',[q(i)*1e-3 rho(i) YFuel(i)-YFuel_st]);
end
fclose(fileID);