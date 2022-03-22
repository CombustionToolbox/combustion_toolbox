function variable = load_struct(filename, variable_name)
    % Load variable from a struct saved in a file
    data = load(filename);
    variable = data.(variable_name);
end