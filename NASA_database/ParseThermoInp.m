function strMaster = ParseThermoInp(reducedDB)
if ~exist('strMaster', 'var')
    if exist('strMaster.mat', 'file') && ~reducedDB
        fprintf('Loading NASA database ... ')
        load('strMaster.mat' , 'strMaster');
    elseif exist('strMaster_reduced.mat', 'file') && reducedDB
        fprintf('Loading Reduced NASA database ... ')
        load('strMaster_reduced.mat' , 'strMaster');
    else
        strMaster = get_strMaster();
        if reducedDB
            strMaster = strMaster_reduced(strMaster);
        end
    end
    fprintf('OK!\n');
else
    fprintf('NASA database already loaded\n')
end
end

function strMaster = get_strMaster()
fid=fopen('thermo.inp'); % loads full database
clc
fprintf('Loading NASA database ... ')
ctLine=0;
ctRefElm=0;
while ctLine<2500
    
    tline = fgetl(fid);
    if ~ischar(tline)
        break
    end
    %     disp(tline)
    if tline(1)=='!'
        continue
    end
    if contains(tline,'thermo')
        tline = fgetl(fid);
        continue
    end
    if contains(tline,'END')
        %         tline = fgetl(fid);
        continue
    end
    ctLine=ctLine+1;
    %% Must have got to a real section
    str.FullName =  sscanf(tline(1:16),'%s');
    str.name = FullName2name(str.FullName);
    str.comments = tline(19:end);
    tline = fgetl(fid);
    str.ctTInt = str2double(tline(1:2));
    str.txRefCode = tline(4:9);
    str.txFormula = tline(11:50);
    str.swtCondensed = str2double(tline(51:52));
    str.mm = str2double(tline(53:65));
    str.Hf0 = str2double(tline(66:80)); % J/mol
    
    if str.ctTInt ==0
        
        tline = fgetl(fid);
        str.tRange = str2double(tline(1:22));
        str.tExponents = str2double(tline(24:63));
        str.Hf298Del0 = str2double(tline(66:end));
    end
    for ctInterval=1:str.ctTInt
        tline = fgetl(fid);
        str.tRange{ctInterval} = str2double(tline(1:22));
        str.tExponents{ctInterval} = str2double(tline(24:63));
        str.Hf298Del0{ctInterval} = str2double(tline(66:end));
        
        tline = fgetl(fid);
        a1 = str2double(tline(1:16));
        a2 = str2double(tline((1:16)+16));
        a3 = str2double(tline((1:16)+32));
        a4 = str2double(tline((1:16)+48));
        a5 = str2double(tline((1:16)+64));
        tline = fgetl(fid);
        a6 = str2double(tline(1:16));
        a7 = str2double(tline((1:16)+16));
        a8 = 0; %str2double(tline((1:16)+32));
        b1 = str2double(tline((1:16)+48));
        b2 = str2double(tline((1:16)+64));
        str.a{ctInterval}=[a1 a2 a3 a4 a5 a6 a7 a8];
        str.b{ctInterval}=[b1 b2];
    end
    
    strMaster.(str.name)=str;
    clear str
end

fclose(fid);
end