function displaysweepresults(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISPLAY SWEEP RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%   varangin = [str,phi,display_species]
%   str = Prop. of state x (phi,species,...)
%   phi = equivalence ratio [-]
%   display_species = display selected species 
% OUTPUT:
%   results on command window
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% help displaysweepresults
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nfrec = 3;
linewidth = 2;
fontsize = 24;
if nargin < 4
    error('Error: Not enough input arguments. Function: displaysweepresults.');
elseif nargin == 5
    display_species = varargin{5};
end
str = varargin{1};
phi = varargin{2};
NameSpecies = varargin{3};
mintol = varargin{4};

if length(phi)>1
% Plot configuration
% color = colours;
f = figure;
set(f,'units','normalized','innerposition',[0.05 0.05 0.9 0.9],...
        'outerposition',[0.05 0.05 0.9 0.9]);
set(axes,'LineWidth',linewidth-0.2,'FontSize',fontsize+2,'BoxStyle','full')
grid on; box on; hold on; axis tight
xlim([min(phi),max(phi)])
ylim([-inf 1])
% ylim([mintol,1])
xlabel('Equivalence Ratio $\phi$','FontSize',fontsize+10,'interpreter','latex');
ylabel('Molar fraction $X_i$','FontSize',fontsize+10,'interpreter','latex');

Nstruct = length(str);
Nspecies = length(str{1}.Xi);
Xi_phi = zeros(length(str{1}.Xi),Nstruct);
if nargin == 4
    for i=1:Nstruct
        Xi_phi(:,i) = str{i}.Xi;
        j = str{i}.Xi>mintol;
        ind = find(j>0);
    %     Xminor(i) = sum(str(i).Xi(~j));
    end
    for i=1:sum(j)
        dl = plot(phi,Xi_phi(ind(i),:),'LineWidth',linewidth);
        if mod(i,nfrec)==0
            loc_label = 'right';
        else
            loc_label = 'left';
        end
        label(dl,NameSpecies{ind(i)},'FontSize',fontsize,'location',loc_label);
    end
elseif nargin == 5
    for i=1:Nstruct
        Xi_phi(:,i) = str{i}.Xi;
    end
    j = [];
    for i = 1:length(NameSpecies)
        if any(strcmpi(NameSpecies{i},display_species))
            j = [j, i];
        end    
    end
    for i=1:length(j)
        dl = plot(phi,Xi_phi(j(i),:),'LineWidth',linewidth);
        if mod(i,nfrec)==0
            loc_label = 'right';
        else
            loc_label = 'left';
        end
        label(dl,NameSpecies{j(i)},'FontSize',fontsize,'location',loc_label);
    end
end
leg =NameSpecies(j);

% plot(phi,Xminor,'color','k','LineWidth',1.2);
set(gca,'yscale','log')
ylim([1e-16 1])
% ylim([-inf 1])
% ylim([0 1])
legend(leg,'FontSize',16,'Location','northeastoutside','interpreter','latex');
% tit = strcat('Molar fraction $X_i$');
% title({tit},'Interpreter','latex','FontSize',16);
% filename2 = strcat(fpath,filename);
% saveas(fig,filename2,'epsc');
else
    error('Funtion displaysweepresults - It is necessary at least 2 cases to display the results.')
end