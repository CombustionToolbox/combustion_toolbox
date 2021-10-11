function tInterval = get_interval(species, T, DB)
    for i = 1:DB.(species).ctTInt
        if (T >= DB.(species).tRange{i}(1)) && (T <= DB.(species).tRange{i}(2))
            break
        end
    end
    tInterval = i;
end