q      = Compute_heatrelease(app.PS.strR,app.PS.strP)*1e-6; % MJ/kg
rho    = Compute_density(app.PS.strR);
Y_Fuel = Compute_YFuel(app.PS.strR,app.PS.strR_Fuel);
W      = Compute_W(app.PS.strR,app.PS.strR_Fuel);
phi = app.PD.phi.Value;
% idx0 = 16;
% q = q(idx0:end);
% rho = rho(idx0:end);
% Y_Fuel = Y_Fuel(idx0:end);

% fid=fopen('CH4.txt','w');
% for i=1:length(q)
%     fprintf(fid,'%f %f %f \n', [q(i), rho(i), Y_Fuel(i)]');
% end
% fclose(fid);