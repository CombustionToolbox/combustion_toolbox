function f = fftfreq(n, varargin)
    % Return the Discrete Fourier Transform sample frequencies
    %
    % The returned float array `f` contains the frequency bin centers in cycles
    % per unit of the sample spacing (with zero at the start). For instance, if
    % the sample spacing is in seconds, then the frequency unit is cycles/second.
    % 
    % Given a window length `n` and a sample spacing `d`:
    % 
    % f = [0, 1, ...,   n/2-1,     -n/2, ..., -1] / (d*n)   if n is even
    % f = [0, 1, ..., (n-1)/2, -(n-1)/2, ..., -1] / (d*n)   if n is odd
    %
    % Args:
    %     n (float): Window length
    %
    % Optional Args:
    %     d (float): Sample spacing (inverse of the sampling rate) [default: 1]
    %
    % Returns:
    %     value (float): Array of length `n//2 + 1` containing the sample frequencies
    
    % Import packages
    import combustiontoolbox.utils.math.floorDiv

    % Definitions
    d = 1;
    n = round(n);

    % Unpack additional inputs
    if (nargin > 1)
        d = varargin{1};
    end
    
    % Calculations
    val = 1.0 / (n * d);
    N = floorDiv(n - 1, 2) + 1;
    p1 = 0:N - 1;
    p2 = -floorDiv(n, 2):-1;
    f = [p1, p2] * val;
end