function DB_master = generate_DB_master(varargin)
    % Generate Database MASTER
    reducedDB = varargin{1,1};
    if nargin > 1
    thermoFile = varargin{1,2};
    else
    thermoFile = 'thermo.inp'; 
    end

    if ~exist('DB_master', 'var')
        if exist('DB_master.mat', 'file') && ~reducedDB
            fprintf('Loading NASA database ... ')
            load('DB_master.mat' , 'DB_master');
        elseif exist('DB_master_reduced.mat', 'file') && reducedDB
            fprintf('Loading Reduced NASA database ... ')
            load('DB_master_reduced.mat' , 'DB_master');
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

    function DB_master = get_DB_master(thermoFile)
    fid=fopen(thermoFile); % loads full database
    clc
    switch thermoFile
        case 'thermo.inp'
            msg = 'Loading NASA database ... ';
        case 'thermo_explo.inp'
            msg = 'Loading NASA database + Explosives database ... ';
        otherwise
            msg = 'Loading an unkown database ... ';
    end

    fprintf(msg)
    ctLine=0;
    ctRefElm=0;
    while ctLine<2500
        
        tline = fgetl(fid);
        if ~ischar(tline)
            break
        end
        if tline(1)=='!'
            continue
        end
        if contains(tline,'thermo')
            tline = fgetl(fid);
            continue
        end
        if contains(tline,'END')
            continue
        end
        ctLine=ctLine+1;
    
        str.FullName =  sscanf(tline(1:16),'%s');
        str.name = FullName2name(str.FullName);
        str.comments = tline(19:end);
        tline = fgetl(fid);
        str.ctTInt = str2double(tline(1:2));
        str.txRefCode = tline(4:9);
        str.txFormula = tline(11:50);
        str.swtCondensed = str2double(tline(51:52));
        str.mm = str2double(tline(53:65));
        str.Hf0 = str2double(tline(66:80));
        
        if str.ctTInt ==0
            tline = fgetl(fid);
            str.tRange = str2num(tline(1:22));
            str.tExponents = str2num(tline(24:63));
            str.Hf298Del0 = str2double(tline(66:end));
        end
        for ctInterval=1:str.ctTInt
            tline = fgetl(fid);
            str.tRange{ctInterval} = str2num(tline(1:22));
            str.tExponents{ctInterval} = str2num(tline(24:63));
            str.Hf298Del0{ctInterval} = str2double(tline(66:end));
            
            tline = fgetl(fid);
            a1 = str2num(tline(1:16));
            a2 = str2num(tline((1:16)+16));
            a3 = str2num(tline((1:16)+32));
            a4 = str2num(tline((1:16)+48));
            a5 = str2num(tline((1:16)+64));
            tline = fgetl(fid);
            a6 = str2num(tline(1:16));
            a7 = str2num(tline((1:16)+16));
            a8 = 0; %str2num(tline((1:16)+32));
            b1 = str2num(tline((1:16)+48));
            b2 = str2num(tline((1:16)+64));
            str.a{ctInterval}=[a1 a2 a3 a4 a5 a6 a7 a8];
            str.b{ctInterval}=[b1 b2];
        end
        
        DB_master.(str.name)=str;
        clear str
    end

    fclose(fid);
end