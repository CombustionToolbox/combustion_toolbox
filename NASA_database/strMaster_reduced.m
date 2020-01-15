function strDB = strMaster_reduced(strDB_0)
NameSpecies = fieldnames(strDB_0);
NSpecies    = numel(NameSpecies);
idx = ones(NSpecies,1);
pattern = {'plus','minus','AL','Ag','F','CL','B','Ca','I','K',...
    'Li','M','D','S','Rb','Pb','V','W','Z','G','T','Cd','Co',...
    'Cr','Cs','Cu','Ni','U','Na','Nb','Hg','CP','HP'};
for i=1:NSpecies
    if contains(NameSpecies(i),pattern)
        idx(i) = 0;
    elseif startsWith(NameSpecies(i),'P')
        idx(i) = 0;
    else
        strDB.(NameSpecies{i}) = strDB_0.(NameSpecies{i});
    end
end

