function DB_master = generate_DB_master(varargin)
    % Generate Mater Database (DB_master) with the thermodynamic data of
    % the chemical species
    %
    % Optional args:
    %     * reducedDB (flag):  Flag indicating reduced database
    %     * thermoFile (file): File with NASA's thermodynamic database
    %
    % Returns:
    %     DB_master (struct):  Database with the thermodynamic data of the chemical species

    reducedDB = varargin{1, 1};

    if nargin > 1
        thermoFile = varargin{1, 2};
    else
        thermoFile = 'thermo.inp';
    end

    if ~exist('DB_master', 'var')

        if exist('DB_master.mat', 'file') && ~reducedDB
            fprintf('Loading NASA database ... ')
            load('DB_master.mat', 'DB_master');
        elseif exist('DB_master_reduced.mat', 'file') && reducedDB
            fprintf('Loading Reduced NASA database ... ')
            load('DB_master_reduced.mat', 'DB_master');
        else
            DB_master = get_DB_master(thermoFile);

            if reducedDB
                DB_master = generate_DB_master_reduced(DB_master);
            end

        end

        fprintf('OK!\n');
    else
        fprintf('NASA database already loaded\n')
    end

end

% SUB-PASS FUNCTIONS
function DB_master = get_DB_master(thermoFile)
    fid = fopen(thermoFile); % loads full database
    clc

    switch thermoFile
        case 'thermo.inp'
            msg = 'Loading NASA database ... ';
        otherwise
            msg = 'Loading an unkown database ... ';
    end

    fprintf(msg)
    line = 0;

    while line < 1e4
        tline = fgetl(fid);

        if ~ischar(tline)
            break
        end

        if isempty(regexp(tline, '\S', 'once'))
            continue
        end

        if tline(1) == '!'
            continue
        end

        if contains(tline, 'thermo')
            tline = fgetl(fid);
            continue
        end

        if contains(tline, 'END')
            continue
        end

        line = line + 1;
        temp.FullName = sscanf(tline(1:16), '%s');
        temp.name = FullName2name(temp.FullName);
        temp.comments = tline(19:end);
        tline = fgetl(fid);
        temp.ctTInt = str2double(tline(1:2));
        temp.txRefCode = tline(4:9);
        temp.txFormula = tline(11:50);
        temp.phase = str2double(tline(51:52));
        temp.mm = str2double(tline(53:65));
        temp.Hf0 = str2double(tline(66:80));

        if temp.ctTInt == 0
            tline = fgetl(fid);
            temp.tRange = str2num(tline(1:22)); %#ok<ST2NM>
            temp.tExponents = str2num(tline(24:63)); %#ok<ST2NM>
            temp.Hf298Del0 = str2double(tline(66:end));
        end

        for ctInterval = 1:temp.ctTInt
            tline = fgetl(fid);
            temp.tRange{ctInterval} = str2num(tline(1:22)); %#ok<ST2NM>
            temp.tExponents{ctInterval} = str2num(tline(24:63)); %#ok<ST2NM>
            temp.Hf298Del0{ctInterval} = str2double(tline(66:end));

            tline = fgetl(fid);
            a1 = str2num(tline(1:16)); %#ok<ST2NM>
            a2 = str2num(tline((1:16) + 16)); %#ok<ST2NM>
            a3 = str2num(tline((1:16) + 32)); %#ok<ST2NM>
            a4 = str2num(tline((1:16) + 48)); %#ok<ST2NM>
            a5 = str2num(tline((1:16) + 64)); %#ok<ST2NM>

            tline = fgetl(fid);
            a6 = str2num(tline(1:16)); %#ok<ST2NM>
            a7 = str2num(tline((1:16) + 16)); %#ok<ST2NM>
            a8 = 0; %str2num(tline((1:16)+32));
            b1 = str2num(tline((1:16) + 48)); %#ok<ST2NM>
            b2 = str2num(tline((1:16) + 64)); %#ok<ST2NM>
            temp.a{ctInterval} = [a1 a2 a3 a4 a5 a6 a7 a8];
            temp.b{ctInterval} = [b1 b2];
        end

        DB_master.(temp.name) = temp;
        clear temp
    end

    fclose(fid);
end
