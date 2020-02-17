function vector = cell2vector(str,field)
    Nstruct = length(str);
    for i=Nstruct:-1:1
        vector(i) = str{i}.(field);
    end
end