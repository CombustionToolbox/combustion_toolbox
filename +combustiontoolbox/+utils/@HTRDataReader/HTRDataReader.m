classdef HTRDataReader < handle
    % The :mat:func:`HTRDataReader` class provides static methods to read 3D fields,
    % coordinates, and general data from HTR solver output files, which are
    % stored in the Hierarchical Data Format (HDF) file format.
    %
    % The Hypersonic Task-based Research (HTR) solver [1] is an open-source
    % exascale-oriented task-based multi-GPU high-order code for hypersonics
    % aerothermodynamics. The HTR solver is available at the GitHub repository [2].
    %
    % References:
    %   [1] Di Renzo, M., Fu, L., and Urzay, J., HTR solver: An open-source
    %       exascale-oriented task-based multi-GPU high-order code for
    %       hypersonics aerothermodynamics, Comput. Phys. Commun, Vol. 255,
    %       2020, p. 107262
    %   [2] https://github.com/stanfordhpccenter/HTR-solver

    methods(Static)

        function value = read(filepath, property)
            % Read general data from the .hdf file
            %
            % Args:
            %     filepath (char): Path to the .hdf file
            %     property (char): Name of the property to read
            %
            % Returns:
            %     value (float): Data array for the specified property
            
            value = h5read(filepath, ['/', property]);
            sz = size(value);
            N_min = min(sz);

            % Transpose if needed
            if sz(1) ~= N_min
                value = value'; 
            end

        end

        function [x, y, z, sz] = read3D(filepath, property)
            % Read a 3D property field from the .hdf file
            %
            % Args:
            %     filepath (char): Path to the .hdf file
            %     property (char): Name of the property to read
            %
            % Returns:
            %     x, y, z (float): 3D arrays with the x, y, z components
            %     sz (float): Size of the 3D array
            
            value = h5read(filepath, ['/', property]);
            sz = size(value);
            value = reshape(value, sz); % Reshape to original dimensions
            sz = sz(2:end); % Extract spatial dimensions
            
            % Extract components
            x = squeeze(value(1, :, :, :));
            y = squeeze(value(2, :, :, :));
            z = squeeze(value(3, :, :, :));
        end

        function [x, y, z, L] = readCoordinates(filepathNodes)
            % Read coordinates and domain length from the .hdf file
            %
            % Args:
            %     filepathNodes (char): Path to the nodes file
            %
            % Returns:
            %     x, y, z (float): 1D arrays of coordinates
            %     L (float): Domain length
            
            % Import packages
            import combustiontoolbox.utils.HTRDataReader;

            [coordinates_x, coordinates_y, coordinates_z, sz_coordinates] = ...
                HTRDataReader.read3D(filepathNodes, 'centerCoordinates');
            
            x = reshape(coordinates_x(:, 1, 1), 1, sz_coordinates(1));
            y = reshape(coordinates_y(1, :, 1), 1, sz_coordinates(2));
            z = reshape(coordinates_z(1, 1, :), 1, sz_coordinates(3));
            
            % Get domain length
            L = max(x);
        end

        function [velocity, sz] = readVelocityField(filepath, property)
            % Read a 3D property field from the .hdf file
            %
            % Args:
            %     filepath (char): Path to the .hdf file
            %     property (char): Name of the property to read
            %
            % Returns:
            %     velocity (VelocityField): Velocity field object
            %     sz (float): Size of the 3D array
            
            % Import packages
            import combustiontoolbox.utils.HTRDataReader;
            import combustiontoolbox.turbulence.VelocityField;

            % Read 3D field
            [u, v, w, sz] = HTRDataReader.read3D(filepath, property);

            % Convert to VelocityField object
            velocity = VelocityField(u, v, w);
        end

    end
    
end
