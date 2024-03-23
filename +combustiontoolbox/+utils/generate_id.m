function id = generate_id(value)
    % Generate a deterministic UUID (Universally Unique Identifier) based
    % on the input character array
    %
    % Args:
    %    value (char): Input character array used to generate the UUID
    %
    % Returns:
    %    id (char): UUID generated from the input character array
    %
    % Example:
    %    id = generate_id('this_is_an_example')
    
    id = java.util.UUID.nameUUIDFromBytes(uint8(value));
end