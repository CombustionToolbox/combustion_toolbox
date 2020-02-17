function display_Tp_phi_2(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISPLAY SWEEP RESULTS (TEMPERATURE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%   varangin = [str1,str2,phi]
%   str1 = Prop. of state x (phi,species,...)
%   str2 = Prop. of state x (phi,species,...)
%   phi = equivalence ratio [-]
% OUTPUT:
%   results on command window
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% help displaysweepresults
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
linewidth = 1.8;
fontsize = 18;
if nargin < 2
    error('Function display_Tp_phi - Not enough input arguments. Function: displaysweepresults.');
end
str1 = varargin{1};
str2 = varargin{2};
phi  = varargin{3};
tit  = varargin{4};

%%% CHECK TIT FOR LATEX
tit = strrep(tit,'%','\%');
%%%

% Plot configuration
% color = colours;
f = figure;
set(f,'units','normalized','innerposition',[0.1 0.1 0.9 0.8],...
    'outerposition',[0.1 0.1 0.9 0.8])
set(axes,'LineWidth',linewidth-0.2,'FontSize',fontsize+2,'BoxStyle','full')
grid on; box on; hold on; axis tight

xlabel('Equivalence Ratio $\phi [-]$','FontSize',fontsize+4,'interpreter','latex');
ylabel('Temperature $T_P [K]$','FontSize',fontsize+4,'interpreter','latex');

Nstruct = length(str1);
for i=Nstruct:-1:1
    Tp1_phi(i) = str1{i}.T;
    Tp2_phi(i) = str2{i}.T;
end
ylim([min(min(Tp1_phi,Tp2_phi)) 1.02*max(max(Tp1_phi,Tp2_phi))])

plot(phi,Tp1_phi,'LineWidth',linewidth);
plot(phi,Tp2_phi,'LineWidth',linewidth);
leg = {'Complete','Incomplete'};
legend(leg,'FontSize',fontsize,'Location','northeast','interpreter','latex');

% tit = strcat('Molar fraction $X_i$');
title({tit},'Interpreter','latex','FontSize',18);
movegui(f,'center')
% filename2 = strcat(fpath,filename);
% saveas(fig,filename2,'epsc');