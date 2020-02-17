function display_Tp_phi(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISPLAY SWEEP RESULTS (TEMPERATURE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%   varangin = [str,phi]
%   str = Prop. of state x (phi,species,...)
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
str = varargin{1};
phi = varargin{2};
tit = varargin{3};
leg = {varargin{4}};

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
% ylim([mintol_display,1])
xlabel('Equivalence Ratio $\phi [-]$','FontSize',fontsize+4,'interpreter','latex');
ylabel('Temperature $T_P [K]$','FontSize',fontsize+4,'interpreter','latex');

Nstruct = length(str);
for i=Nstruct:-1:1
    Tp_phi(i) = str{i}.T;
end
plot(phi,Tp_phi,'LineWidth',linewidth);
% tit = strcat('Molar fraction $X_i$');
legend(leg,'FontSize',fontsize,'Location','northeast','interpreter','latex');
title({tit},'Interpreter','latex','FontSize',16);
movegui(f,'center')
% filename2 = strcat(fpath,filename);
% saveas(fig,filename2,'epsc');