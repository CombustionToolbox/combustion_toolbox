% load DET_H2_fixed2;

[rho, ~, W] = plot_density_YFuel(app.PS.strR,app.PS.strR_Fuel);
[q, YFuel, H, YFuel_st] = plot_heatrelease_YFuel(app.PS.strR,app.PS.strP,app.PS.strR_Fuel);

points = find_common(app.PD.phi.Value, W, H);
fileID = fopen(strcat('C',sprintf('%d',app.PS.strR{1}.x),'H',sprintf('%d',app.PS.strR{1}.y),'.txt'),'w');
for i=1:length(rho)
    fprintf(fileID,'%f   %f   %f\n',[q(i)*1e-3 rho(i) YFuel(i)-YFuel_st]);
end
fclose(fileID);

fileID = fopen(strcat('HW_C',sprintf('%d',app.PS.strR{1}.x),'H',sprintf('%d',app.PS.strR{1}.y),'.txt'),'w');
% for i=1:25:length(H)
%     fprintf(fileID,'%f   %f   %f\n',[H(i) W(i) H(i)/W(i)]);
% end

fprintf(fileID,'%f   %f   %f\n',[H(21) W(21) H(21)/W(21)]);
fprintf(fileID,'%f   %f   %f\n',[H(61) W(61) H(61)/W(61)]);
fprintf(fileID,'%f   %f   %f\n',[H(101) W(101) H(101)/W(101)]);
fprintf(fileID,'%f   %f   %f\n',[H(141) W(141) H(141)/W(141)]);
fprintf(fileID,'%f   %f   %f\n',[H(181) W(181) H(181)/W(181)]);
fprintf(fileID,'%f   %f   %f\n',[H(221) W(221) H(221)/W(221)]);
fprintf(fileID,'%f   %f   %f\n',[H(261) W(261) H(261)/W(261)]);
fprintf(fileID,'%f   %f   %f\n',[H(301) W(301) H(301)/W(301)]);
fprintf(fileID,'%f   %f   %f\n',[H(341) W(341) H(341)/W(341)]);

fclose(fileID);