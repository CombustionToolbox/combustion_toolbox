function id = generateID(value)
    % Generate a deterministic numeric identifier from input data.
    %
    % Args:
    %    value (char): Input character array used to generate the ID
    %
    % Returns:
    %    id (double): Deterministic 32-bit numeric identifier
    %
    % Example:
    %    id = generateID('this_is_an_example')

    % Compute 128-bit MD5 hash from input string
    md = java.security.MessageDigest.getInstance('MD5');
    md.update( uint8(value) );

    % Convert the 16-byte MD5 hash into two uint32 integers
    hash = typecast(md.digest, 'uint32');

    % Collapse the 128-bit digest into a single 32-bit numeric identifier
    id = double( bitxor( bitxor( hash(1), hash(2) ), bitxor( hash(3), hash(4) ) ) );
end