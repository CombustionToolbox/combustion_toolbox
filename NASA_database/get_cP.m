function cP = get_cP(species, T, DB)
    if DB.(species).ctTInt > 0
        a = DB.(species).a;
        b = DB.(species).b;
        tRange       = DB.(species).tRange;
        tExponents   = DB.(species).tExponents;
        ctTInt       = DB.(species).ctTInt;
        R0           = 8.3144598;
        if (T < tRange{1}(1)) || (T > tRange{ctTInt}(2))
            cP = species_cP(species, T, DB);
            return
        else
            for i = 1:ctTInt
                if (T >= tRange{i}(1)) && (T <= tRange{i}(2))
                    tInterval = i;
                end
            end
            cP = R0 * sum(a{tInterval} .* T.^tExponents{tInterval});
        end
    end
end
