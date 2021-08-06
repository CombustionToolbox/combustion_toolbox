function cP = get_cP(species, T, strDB)
    if strDB(species).ctTInt > 0
        a = strDB.(Species).a;
        b = strDB.(Species).b;
        tRange       = strDB.(Species).tRange;
        tExponents   = strDB.(Species).tExponents;
        
        if (T < tRange{1}(1)) || (T > tRange{ctTInt}(2))
            cP = species_cP(species, T, strDB);
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
