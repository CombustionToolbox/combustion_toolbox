function tInterval = compute_interval_NASA(species, T, DB, tRange, ctTInt)
    % Compute interval NASA polynomials
    
    tInterval = [];
    if DB.(species).ctTInt > 0
        for j = 1:ctTInt
            if (T >= tRange{j}(1)) && (T <= tRange{j}(2))
                tInterval = j;
            end
        end
        if isempty(tInterval)
            if  T > tRange{j}(2)
                tInterval = ctTInt;
            elseif T < tRange{j}(2)
                tInterval = 1;
            end
        end
    end
end