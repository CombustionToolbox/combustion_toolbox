function q = Compute_heatrelease(strR,strP)
if length(strR)>1
    for i=length(strR):-1:1
        mR(i)   = strR{i}.mi;
        dhR(i)  = strR{i}.DhT*1e3/mR(i); %J/kg == m^2/s^2
        uR(i)   = strR{i}.u;
        mP(i)   = strR{i}.mi;
        dhP(i)  = strP{i}.DhT*1e3/mP(i); %J/kg == m^2/s^2
        uP(i)   = strP{i}.u;
        % aR    = strR.sound;
        % gamma = strP.gamma;
    end
else
    mR   = strR.mi;
    dhR  = strR.DhT*1e3/mR; %J/kg == m^2/s^2
    uR   = strR.u;
    mP   = strR.mi;
    dhP  = strP.DhT*1e3/mP; %J/kg == m^2/s^2
    uP   = strP.u;
    % aR    = strR.sound;
    % gamma = strP.gamma;
end
q = (dhP-dhR+0.5*(uP.^2-uR.^2));
% q = q.*(gamma.^2-1)./(2*aR.^2);