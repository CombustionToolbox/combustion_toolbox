function rho = Compute_density(mix)
    % Get density of the set of mixtures
    %
    % Args:
    %     mix (struct): Properties of the mixture/s
    %
    % Returns:
    %     rho (float):  Vector with the densities of all the mixtures

    for i=length(mix):-1:1
        rho(i) = mix{i}.rho;
    end
end