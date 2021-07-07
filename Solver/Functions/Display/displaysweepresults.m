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

if nargin < 5
    error('Error: Not enough input arguments. Function: displaysweepresults.');
end
app = varargin{1};
str = varargin{2};
x = varargin{3};
NameSpecies = varargin{4};
mintol = varargin{5};
config = varargin{6}; 
NE = length(NameSpecies);
if nargin == 7
    display_species = varargin{7}; 
end
if length(x)>1
    % Plot configuration
    % color = colours;
    if ~app.Misc.FLAG_GUI
        f = figure;
        set(f,'units','normalized','innerposition',[0.05 0.05 0.9 0.9],...
            'outerposition',[0.05 0.05 0.9 0.9]);
        set(axes,'LineWidth',config.linewidth,'FontSize',config.fontsize-2,'BoxStyle','full')
    else
        f = app.UIAxes;
        cla(f);
    end
    grid(f, 'off'); box(f, 'off'); hold(f, 'on'); 
    
    xlabel(f, config.labelx,'FontSize',config.fontsize,'interpreter','latex');
    ylabel(f, config.labely,'FontSize',config.fontsize,'interpreter','latex');
    
    Nstruct = length(str);
%     Nspecies = length(str{1}.Xi);
    Xi_phi = zeros(length(str{1}.Xi),Nstruct);
    if nargin == 6
%         for i=1:Nstruct
%             Xi_phi(:,i) = str{i}.Xi;
%         end
%         j = any(Xi_phi'>mintol);
%         
%         NE = sum(j);
        all_ind = [];
        for i=1:Nstruct
            Xi_phi(:,i) = str{i}.Xi;
            j = str{i}.Xi>mintol;
            ind = find(j>0);
            all_ind = [all_ind; ind];
            %     Xminor(i) = sum(str(i).Xi(~j));
        end
        all_ind = unique(all_ind);
        NE = numel(all_ind);  
        colorbw = brewermap(NE,'Spectral');
        for i=1:numel(all_ind)
            dl = plot(f, x, Xi_phi(all_ind(i),:),'LineWidth',config.linewidth,'color',colorbw(i,:));
%             if mod(i,nfrec)==0
%                 loc_label = 'right';
%             else
%                 loc_label = 'left';
%             end
%             label(dl,NameSpecies{ind(i)},'FontSize',fontsize,'location',loc_label);
        end
    elseif nargin == 7
        for i=1:Nstruct
            Xi_phi(:,i) = str{i}.Xi;
        end
        all_ind = [];
        for i = 1:length(NameSpecies)
            if any(strcmpi(NameSpecies{i},display_species))
                all_ind = [all_ind, i];
            end
        end
        k = any(Xi_phi(all_ind,:)'>mintol);
        all_ind = all_ind(k);
        
        NE = numel(all_ind);  
        colorbw = brewermap(NE,'Spectral');
        for i=1:numel(all_ind)
            dl = plot(f, x,Xi_phi(all_ind(i),:),'LineWidth',config.linewidth,'color',colorbw(i,:));
%             if mod(i,nfrec)==0
%                 loc_label = 'right';
%             else
%                 loc_label = 'left';
%             end
%             label(dl,NameSpecies{j(i)},'FontSize',fontsize,'location',loc_label);
        end
    end
    leg =NameSpecies(all_ind);
    
    % plot(phi,Xminor,'color','k','LineWidth',1.2);
    set(f,'yscale','log')
    
%     xlim(f, [min(x),max(x)])
    
    ymin = 10^floor(log(abs(min(Xi_phi(Xi_phi>0))))/log(10));
    if ymin > mintol
        ylim(f, [ymin,1])
    else
        ylim(f, [mintol,1])
    end

    legend(f, leg,'FontSize',config.fontsize-6,'Location','northeastoutside','interpreter','latex');
    % tit = strcat('Molar fraction $X_i$');
    % title({tit},'Interpreter','latex','FontSize',16);
    % filename2 = strcat(fpath,filename);
    % saveas(fig,filename2,'epsc');
%     movegui(f,'center')
else
    error('Funtion displaysweepresults - It is necessary at least 2 cases to display the results.')
end