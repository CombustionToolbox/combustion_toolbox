function WorH = plot_WorH_YFuel(str,x,y,Wbool)
linewidth = 2;
fontsize = 24;

dy0 = ones(1,length(y));

h = x(2:end)-x(1:end-1);
h(end+1) = h(end);
hmax = max(h);

for i=length(h):-1:1
    mu(i) = h(i)/hmax;
end
for i=length(y)-1:-1:2
    dy0(i) = 1/y(i)*(y(i+1)-y(i-1))/((mu(i+1)+mu(i))*hmax);
end
dy0(1)   = 1/y(1)*(y(2)-y(1))/(x(2)-x(1));
dy0(end) = 1/y(end)*(y(end)-y(end-1))/(x(end)-x(end-1));


if ~isfield(str{1},'z')
    str{1}.z = 0;
end
s = str{1}.x + str{1}.y/4 - str{1}.z/2;

WA = 28.850334000000007;
WF = str{1}.x*12.0107 + str{1}.y*1.0079 + str{1}.z*15.9994;
for i=length(str):-1:1
    mo(i) = s*(2*15.9994 + 79/21*2*14.0067)/str{i}.phi;
    rho_Fuel(i) = str{i}.rho;
end
Nf = 1;
mu = Nf*WF;
m  = mu + mo;
Yu = mu./m;

rhoA = 100000/((8314/WA)*300);
rhou = 1./(1-(1-WA/WF).*Yu)*rhoA;
W = (1-WA/WF)./(1-(1-WA/WF).*Yu);

f = figure;
set(f,'units','normalized','innerposition',[0.05 0.05 0.9 0.9],...
        'outerposition',[0.05 0.05 0.9 0.9]);
set(axes,'LineWidth',linewidth-0.2,'FontSize',fontsize+2,'BoxStyle','full')
movegui(f,'center')
grid on; box on; hold on; axis tight
xlabel('$Y_{Fuel}\ [-]$','FontSize',fontsize+10,'interpreter','latex');
if Wbool
    ylabel('$W$','FontSize',fontsize+10,'interpreter','latex');
else
    ylabel('$H$','FontSize',fontsize+10,'interpreter','latex');
end
plot(x,dy0,'LineWidth',linewidth)

if Wbool
    plot(x,W,'LineWidth',linewidth)
end

WorH = dy0;
