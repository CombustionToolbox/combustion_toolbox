function [str1,str2,guess] = SolveProblemDET_OVERDRIVEN(app, strR, p1, T1, guess, overdriven)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE CHAPMAN-JOUGET STATE (CJ UPPER STATE - STRONG DETONATION)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% help SolveProblemDET
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GLOBAL CONFIGURATION
counter = 1;
numsteps = 5;
max = guess(4);
min = guess(5);
Goodness = 0.99999;
g = 0;

% guess(1) = 2000;
% guess(2) = 2000;
% str1 = strR;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while  ((counter <= 3) || (g <= Goodness))
    step = (max-min)/numsteps;
    i = 1;
    w1 = zeros(1,numsteps+1);
    rr = zeros(1,numsteps+1);
    for x = min:step:max
        [str1, str2, guess_shock] = SolveProblemDET_inloop(app, strR, p1, T1, x, guess);
        w1(i) = str1.u;
        rr(i) = str2.rho/strR.rho;
        
%         guess(1) = str2.T;
%         guess(2) = str1.u;
        
        i = i + 1;
    end
    
%     F = fitoptions('method','NonlinearLeastSquares','StartPoint',[1,1,1]);
%     ftpe = fittype('a*x^2+b*x+c','coeff',{'a','b','c'});
%     [curve,goodness,~]=fit(transpose(rr),transpose(w1),ftpe,F);
%     
%     a = curve.a;
%     b = curve.b;
%     
%     g = goodness.rsquare;
%     dnew = -b/(2*a);

    F = fitoptions('method','NonlinearLeastSquares','StartPoint',[1,1,1,1]);
    ftpe = fittype('a0*x^3+a*x^2+b*x+c','coeff',{'a0','a','b','c'});
    [curve,goodness,~]=fit(transpose(rr),transpose(w1),ftpe,F);
    
    a0 = curve.a0;
    a = curve.a;
    b = curve.b;
    
    g = goodness.rsquare;
    dnew = (-a+sqrt(a^2-3*a0*b))/(3*a0);
    if ~isreal(dnew)|| dnew<0
        dnew = (-a-sqrt(a^2-3*a0*b))/(3*a0);
    end
    
%     min = dnew - dnew*0.0005;
%     max = dnew + dnew*0.0005;
    if counter <= 3
        min = dnew - dnew*0.001;
        max = dnew + dnew*0.001;
    elseif counter <= 15
        min = dnew - dnew*0.0005;
        max = dnew + dnew*0.0005;
    else
        min = dnew - dnew*0.0001;
        max = dnew + dnew*0.0001;
        % TOLERANCES
%         TN.ERRFT = 1e-4;
%         TN.ERRFU = 1e-4;
        Goodness = 0.9999;
    end
    counter = counter + 1;
end

[~,cj_speed] = predint(curve,dnew,0.95,'obs','on');
% [h,cj_speed] = predint(curve,dnew,0.95,'functional','off');

% U_cj = cj_speed + .001;  %Make sure Hugoniot and Rayleigh intersect
U_cj = cj_speed;  %Make sure Hugoniot and Rayleigh intersect
U_cj = U_cj*overdriven;
clear str1 str2
for i=length(overdriven):-1:1
% guess_shock(3) = dnew;
% TN.guess_shock = guess_shock;
% [str1,str2] = SolveProblemSHOCK_I(strR,phi,p1,T1,U_cj,E,S,C,M,PD,TN,strThProp);
[str1.overdriven{i},str2.overdriven{i}] = SolveProblemSHOCK_Jouget_I(app, strR, p1, T1, U_cj(i));
if i==1
    % guess(1) = str2.T;
    % guess(2) = U_cj;
    guess(4) = min;
    guess(5) = max;
end
end
%%% PLOT CURVE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% linewidth = 1.6;
% fontsize = 24;
% color = colours;
% f = figure;
% x = linspace(dnew-0.01,dnew+0.01);
% 
% set(f,'units','normalized','innerposition',[0.185 0.185 0.55 0.8],...
%         'outerposition',[0.185 0.185 0.75 0.8]);
% set(axes,'LineWidth',linewidth-0.2,'FontSize',fontsize+2,'BoxStyle','full')
%     hold on; grid on; box on; axis tight
% xlabel('$\rho_1/\rho_2$','FontSize',fontsize+4,'interpreter','latex');
% ylabel('Velocidad, $w [m/s]$','FontSize',fontsize+4,'interpreter','latex');
% 
% plot(rr,w1,'d','MarkerFaceColor',color(2,:),'MarkerEdgeColor','black',...
%                 'MarkerSize',fontsize-12,'color',color(2,:));
% plot(x,curve(x),'-','LineWidth',linewidth,'color',color(2,:))
% plot(dnew,U_cj,'s','MarkerFaceColor',color(1,:),'MarkerEdgeColor','black',...
%                 'MarkerSize',fontsize-12,'color',color(1,:));