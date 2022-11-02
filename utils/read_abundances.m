function [abundances, elements] = read_abundances(filename)
    % Read solar abundances file
    % Format: [number element, element, abundance, name, molar mass (g/mol)]
    %
    % Args:
    %    filename (file): Filename with the data
    %
    % Returns:
    %    Tuple containing:
    %
    %    * abundances (float): Vector with the logarithmic base 10 solar abundances
    %    * elements (cell): List with the given elements

    % Open file
    fileID = fopen(filename, 'r');
    % Define the format of the data to read and the shape of the output array
    formatSpec = '%d %s %f %s %f';
    % Read the file data, filling output array, A, in column order. fscanf
    % reuses the format, formatSpec, throughout the file.
    A = textscan(fileID, formatSpec);
    % Close file
    fclose(fileID);
    % Get inputs
    abundances = A{1, 3};
    elements = A{1, 2};
end
