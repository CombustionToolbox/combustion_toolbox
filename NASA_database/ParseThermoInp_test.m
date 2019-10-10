function strMaster = ParseThermoInp_test()
fid=fopen('thermo.inp'); % loads full database
clc
fprintf('Loading NASA database ... ')
ctLine=0;
ctRefElm=0;
while ctLine<2500
    
    tline = fgetl(fid);
    if ~ischar(tline), break, end
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
    str.ctTInt = str2num(tline(1:2));
    str.txRefCode = tline(4:9);
    str.txFormula = tline(11:50);
    str.swtCondensed = str2num(tline(51:52));
    str.mm = str2num(tline(53:65));
    str.Hf0 = str2num(tline(66:80)); % J/mol
    
    if str.ctTInt ==0
        
        tline = fgetl(fid);
        str.tRange = str2num(tline(1:22));       
        str.tExponents = str2num(tline(24:63));
        str.Hf298Del0 = str2num(tline(66:end));
    end
    for ctInterval=1:str.ctTInt
        tline = fgetl(fid);
        str.tRange{ctInterval} = str2num(tline(1:22));
        str.tExponents{ctInterval} = str2num(tline(24:63));
        str.Hf298Del0{ctInterval} = str2num(tline(66:end));
        
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
    
%    Uncomment to list reactants and their reference temperatures
%
%    if ~iscell(str.tRange)
%        if str.tRange(1) < 93 % The U.S. National Institute of Standards and Technology has chosen to consider the field of cryogenics as that involving temperatures below ?180 °C (93 K; ?292 °F). (Wikipedia)
%            disp(['> ',str.FullName,' [',num2str(str.tRange(1)),' - ',num2str(str.tRange(2)),'] Cryogenic reactant'])
%        elseif str.tRange(1) < 298.15
%            disp(['> ',str.FullName,' [',num2str(str.tRange(1)),' - ',num2str(str.tRange(2)),'] Reactant at low T'])
%        else
%            disp(['> ',str.FullName,' [',num2str(str.tRange(1)),' - ',num2str(str.tRange(2)),'] Reactant at room T'])
%        end
%    end
%     
%    Uncomment to list reference form of elements and their respective
%    temperature ranges
%
%    if ~isempty(strfind(str.comments,'Ref-Elm.'))
%        ctRefElm = ctRefElm +1;
%        disp(['''',name_with_parenthesis(str.name),' [',num2str(str.tRange{1}(1)),'-',num2str(str.tRange{str.ctTInt}(2)),']'';'])
%    end
%    
%    disp(tline)
%    disp([str.name ' : ' str.FullName ' : ' str.comments] )
%
%    Uncomment to list all species in NASA database (2002)
%
%    fprintf('%18s ',str.FullName)
%    if mod(ctLine,5)==0
%        fprintf('\n')
%    end
    
    strMaster.(str.name)=str;
    clear str
end

fprintf('OK!\n');

fclose(fid);
