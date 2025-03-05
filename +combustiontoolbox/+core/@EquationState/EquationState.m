classdef (Abstract) EquationState < handle
    % The :mat:func:`EquationState` class is an abstract class that defines
    % the interface for computing the pressure and molar volume of a
    % mixture using a specified equation of state.
    %
    % Subclasses must implement the following abstract methods:
    %     * getPressure(temperature, molarVolume, varargin)
    %     * getVolume(temperature, pressure, varargin)
    %
    % See also: :mat:func:`EquationStateIdealGas`, :mat:func:`Mixture`
    
    methods (Abstract)
        % Compute pressure [Pa] given the temperature and molar volume.
        pressure = getPressure(obj, temperature, molarVolume, varargin)
        
        % Compute molar volume [m3/mol] given the temperature and pressure.
        molarVolume = getVolume(obj, temperature, pressure, varargin)
    end
    
end