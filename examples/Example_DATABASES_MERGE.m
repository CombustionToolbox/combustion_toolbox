% -------------------------------------------------------------------------
% EXAMPLE: DATABASES MERGE
%
% Example of merging two databases. In this case, a NASA database and a
% Burcat database using NASA9 format, which has been previously converted
% using the static method thermoMillennium2thermoNASA9 from the Burcat
% class
%
% @author: Alberto Cuadra Lara
%                 
% Last update October 14 2025
% -------------------------------------------------------------------------

% Import packages
import combustiontoolbox.databases.*

% Get Nasa database
DB_NASA = NasaDatabase('thermoFile', 'thermo_NASA.inp');

% Generate Burcat database using NASA9 format
DB_BURCAT = NasaDatabase('thermoFile', 'thermo_millennium_2_thermoNASA9.inp');

% Merge databases
DB = DB_NASA + DB_BURCAT;