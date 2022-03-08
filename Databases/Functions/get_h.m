function h = get_h(species, T, DB)
    if DB.(species).ctTInt > 0
        R0 = 8.3144598;
        a = DB.(species).a;
        b = DB.(species).b;
        ctTInt       = DB.(species).ctTInt;
        tRange       = DB.(species).tRange;
        tExponents   = DB.(species).tExponents;
        
        if (T < tRange{1}(1)) || (T > tRange{ctTInt}(2))
            h = species_DhT(species, T, DB);
            return
        else
            for i = 1:ctTInt
                if (T >= tRange{i}(1)) && (T <= tRange{i}(2))
                    tInterval = i;
                end
            end
            h = R0 * T * (sum(a{tInterval} .* T.^tExponents{tInterval} .* [-1   log(T) 1      1/2 1/3 1/4 1/5 0]) + b{tInterval}(1)/T);
            h = h*1e-3;
        end
    end
end