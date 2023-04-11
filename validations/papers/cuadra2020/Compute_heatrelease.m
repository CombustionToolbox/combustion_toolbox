function [q, Q] = Compute_heatrelease(mix1, mix2)
    % Compute heat release [J/kg] of the chemical transformation of the mixture 1 to mixture 2
    %
    % Args:
    %     mix1 (struct): Properties of the initial mixture
    %     mix2 (struct): Properties of the final mixture
    %
    % Returns:
    %     q (float): heat release [J/kg] == [m^2/s^2]
    %     Q (float): dimensionless heat release

    if length(mix1)>1
        for i=length(mix1):-1:1
            mR(i)  = mix1{i}.mi;
            dhR(i) = mix1{i}.DhT*1e3 / mR(i); %J/kg == m^2/s^2
            uR(i)  = mix1{i}.u;
            mP(i)  = mix2{i}.mi;
            dhP(i) = mix2{i}.DhT*1e3 / mP(i); %J/kg == m^2/s^2
            uP(i)  = mix2{i}.v_shock;
            aR    = mix1.sound;
            gamma = mix2.gamma;
        end
    else
        mR  = mix1.mi;
        dhR = mix1.DhT*1e3/mR; %J/kg == m^2/s^2
        uR  = mix1.u;
        mP  = mix2.mi;
        dhP = mix2.DhT*1e3/mP; %J/kg == m^2/s^2
        uP  = mix2.v_shock;
        aR  = mix1.sound;
        gamma = mix2.gamma;
    end
    q = (dhP - dhR + 0.5*(uP.^2 - uR.^2));
    Q = q.*(gamma.^2 - 1) ./ (2 * aR.^2);
end