function points = find_common(x, y1, y2)

linewidth = 1.2;
fontsize = 12;
f = figure;
set(f,'units','normalized','innerposition',[0.1 0.1 0.7 0.9],...
        'outerposition',[0.1 0.1 0.7 0.9]);
movegui(f,'center')
set(axes,'LineWidth',linewidth-0.2,'FontSize',fontsize+2,'BoxStyle','full')
grid on; box on; hold on
xlim([min(x),max(x)])
xlabel('Equivalence Ratio, $\phi$','FontSize',fontsize+10,'interpreter','latex');
ylabel('$W$ and $H$','FontSize',fontsize+10,'interpreter','latex');
% xmin = 0.1;
% xmax = 3.5;
% xprime = linspace(xmin, xmax, 200);

% p2 = polyfit(x, y2, 2);
% 
% f2 = polyval(p2, xprime);
% 
% y2prime = log(y2); % take natural logarithm of y data values
% p2prime = polyfit(x,y2prime,1);
% aprime = p2prime(2);
% bprime = p2prime(1);
% 
% for i = 1:length(y1)
% %     f = [p1(1); p1(2); p1(3); p1(4); p1(5)-y2(i)]; % Polynomial == 0 coifficients
%     for z = 1:2
%         if z 
%             p1 = polyfit(x, abs(y1), 2);
%         else
%             p1 = polyfit(x, -abs(y1), 2);
%         end
%         f1 = polyval(p1, xprime);
%         f = [p1(1); p1(2); p1(3)- y2(i)]; % Polynomial == 0 coifficients
%         r1 = roots(f);
%         if ~isempty(r1)
%             k = 1;
%             points = [];
%             for j=length(r1):-1:1
%                 if isreal(r1(j)) && r1(j) >= xmin && r1(j) <= xmax 
%                     points(1, k) = r1(j);
%                     points(2, k) = polyval(p1, r1(j));
%                     k = k + 1; 
%                 end
%             end
%         end
%     end
% end  
l(1) = plot(x, y1, '-b','LineWidth',linewidth);
l(2) = plot(x, -y1, '-r','LineWidth',linewidth);
l(3) = plot(x, y2, '-k','LineWidth',linewidth);
% plot(xprime, f1, '-r')
% plot(xprime, exp(aprime)*exp(xprime*bprime), '-g')

p = InterX([x; y1],[x; y2]);
if ~isempty(p)
    points1 = p;
    l(4) = plot(p(1,:), p(2,:), 'bo', 'markersize', 8);
end
p = InterX([x; -y1],[x; y2]);
if ~isempty(p)
    points2 = p;
    l(5) = plot(p(1,:), p(2,:), 'ro', 'markersize', 8);
end
leg = {'$W$', '-$W$', '$H$', '$W = H$', '$-W = H$'};
legend(l, leg,'Location', 'northeast', 'interpreter', 'latex')

points = [points1, points2];
