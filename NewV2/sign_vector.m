function y = sign_vector(x)
    % sign_vector - Computes the sign of each element in a vector or array.
    %
    % Syntax:
    %   y = sign_vector(x)
    %
    % Input:
    %   x - A numeric vector or array.
    %
    % Output:
    %   y - An array of the same size as x, with each element being:
    %       1  if the corresponding element in x > 0,
    %       -1 if the corresponding element in x < 0,
    %       0  if the corresponding element in x == 0.

    % Validate input
    if ~isnumeric(x)
        error('Input must be a numeric vector or array.');
    end
    
    % Compute the sign for each element
    y = zeros(size(x)); % Initialize the output array
    y(x > 0) = 1;       % Positive values
    y(x < 0) = -1;      % Negative values
    % Zero values remain zero (already initialized as 0)
end
