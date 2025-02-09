function [cosTheta] = getAlignment(field1, field2)
    % Computes the local cosine of the angle (alignment) between two vector fields
    %
    % Args:
    %     field1 (float): 3D vector field
    %     field2 (float): 3D vector field
    %
    % Returns:
    %     cosTheta (float): 3D array representing the pointwise cosine of the angle between two vectors

    % Compute pointwise dot product
    dotProduct = field1 .* field2;

    % Compute magnitudes (norms) of the two fields
    magnitude1 = norm(field1);
    magnitude2 = norm(field2);

    % Compute cosine of angles between gradients
    cosTheta = dotProduct ./ (magnitude1 .* magnitude2);

    % Handle numerical stability (e.g., division by zero)
    cosTheta(~isfinite(cosTheta)) = NaN; % Set invalid values to NaN for clarity
end