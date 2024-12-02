function [cosTheta] = getAlignment(field1, field2)
    % Computes the local cosine of the angle (alignment) between two scalar fields
    %
    % Args:
    %     field1 (float): 3D scalar field
    %     field2 (float): 3D scalar field
    %
    % Returns:
    %     cosTheta (float): 3D array representing the pointwise cosine of the angle between two fields

    % Compute pointwise dot product
    dotProduct = field1 .* field2;

    % Compute pointwise magnitudes
    magnitude1 = abs(field1);
    magnitude2 = abs(field2);

    % Compute cosine of angles between gradients
    cosTheta = dotProduct ./ (magnitude1 .* magnitude2);

    % Handle numerical stability (e.g., division by zero)
    cosTheta(~isfinite(cosTheta)) = NaN; % Set invalid values to NaN for clarity
end