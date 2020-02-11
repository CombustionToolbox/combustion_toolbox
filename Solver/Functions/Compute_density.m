function rho = Compute_density(str)
for i=length(str):-1:1
    rho(i) = str{i}.rho;
end