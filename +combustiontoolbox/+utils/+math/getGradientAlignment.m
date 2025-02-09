function [cosTheta] = getGradientAlignment(field1, field2, x, y, z)
    % Compute the alignment between two gradient fields
    %
    % Args:
    %     field1 (float): 3D scalar field for which the gradient is computed
    %     field2 (float): 3D scalar field for which the gradient is computed
    %     x (float): 3D array with the x-coordinate of the grid
    %     y (float): 3D array with the y-coordinate of the grid
    %     z (float): 3D array with the z-coordinate of the grid
    %
    % Returns:
    %     cosTheta (float): 3D array with the cosine of the angle between gradients

    % Import packages
    import combustiontoolbox.utils.math.gradientPeriodic;

    % Get gradients
    [df1_dx, df1_dy, df1_dz] = gradientPeriodic(field1, x, y, z);
    [df2_dx, df2_dy, df2_dz] = gradientPeriodic(field2, x, y, z);

    % Compute dot product of gradients
    dotProduct = df1_dx .* df2_dx + df1_dy .* df2_dy + df1_dz .* df2_dz;

    % Compute magnitudes of gradients
    magnitude1 = sqrt(df1_dx.^2 + df1_dy.^2 + df1_dz.^2);
    magnitude2 = sqrt(df2_dx.^2 + df2_dy.^2 + df2_dz.^2);

    % Compute cosine of angles between gradients
    cosTheta = dotProduct ./ (magnitude1 .* magnitude2);

    % Handle numerical stability (e.g., division by zero)
    cosTheta(~isfinite(cosTheta)) = NaN; % Set invalid values to NaN for clarity
end