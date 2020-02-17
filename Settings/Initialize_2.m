function [E,S,M,C,Misc,PT] = Initialize_2(E,S,M,C,Misc)
%% Incomplete combustion
% First we eliminate from the minor species list those considered major
% species in case the user has included any
n_pass = [];
for n = 1:M.Lminor
    if ~any(strcmp(M.minor_products{n},S.List_fixed_Species))
        n_pass = [n_pass, n];
    end
end
M.minor_products = M.minor_products(n_pass);
S.List_Compute_Species = [S.List_fixed_Species,M.minor_products];
%%
if Misc.FLAG_FIRST
    M.L_minor = length(M.minor_products);
    if M.L_minor > 0
        for n=M.L_minor:-1:1
            % Properties of other minor species under consideration, which can
            % be written in the generic form C_alpha H_beta O_gamma N_omega
            % Find index minor species
            M.idx_minor(n) = find(strcmp(S.NameSpecies,M.minor_products{n}));
            Misc.FLAG_FIRST = false;
        end
        C.alpha = C.A0.Value(M.idx_minor,E.ind_C)';
        C.beta  = C.A0.Value(M.idx_minor,E.ind_H)';
        C.gamma = C.A0.Value(M.idx_minor,E.ind_O)';
        C.omega = C.A0.Value(M.idx_minor,E.ind_N)';
        
        S.idx_all = [S.idx_fixed,M.idx_minor];
    else
        S.idx_all = S.idx_fixed;
    end
end
%% CH4 major specie
if any(contains(M.minor_products,'CH4'))
    M.major_CH4 = true;
    M.idx_m_CH4 = find_idx({'CH4'},M.minor_products);
%     M.idx_m_C2H2= find_idx({'C2H2_acetylene'},M.minor_products);
%     M.idx_m_C6H6= find_idx({'C6H6'},M.minor_products);
    
    M.idx_m_CH3 = find_idx({'CH3'},M.minor_products);
    M.idx_m_H = find_idx({'H'},M.minor_products);
    M.idx_m_CH = find_idx({'CH'},M.minor_products);
%     M.idx_m_C = find_idx({'C'},M.minor_products);
else
    M.major_CH4 = false;
end
%% OH major specie
if any(contains(M.minor_products,'OH'))
    M.major_OH = true;
    M.idx_m_OH = find_idx({'OH'},M.minor_products);
else
    M.major_OH = false;
end
%% C minor specie
% if any(contains(M.minor_products,'C'))
%     M.minor_C = true;
%     M.idx_m_C = find_idx({'C'},M.minor_products);
% else
%     M.minor_C = false;
% end
%% Ask problem type
% prompt = {'Enter problem type:'};
% dlgtitle = 'Input';
% dims = [1 35];
% definput = {'HP'};
% PT = inputdlg(prompt,dlgtitle,dims,definput);
% PT = PT{1};

fn = {'TP','HP','SP','TV','EV','SV','SHOCK_I','SHOCK_R','DET','DET_OVERDRIVEN'};
[indx,tf] = listdlg('PromptString','Select a problem:',...
                           'SelectionMode','single',...
                           'ListString',fn,'ListSize',[150,120]);
if tf
    PT = fn{indx};
else
    error('Problem type not selected.')
end
% output = fopen('output.txt','wt');
% f = waitbar(0,'Please wait...','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
% setappdata(f,'canceling',0);