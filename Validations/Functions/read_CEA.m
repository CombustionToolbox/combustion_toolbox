% READ SCRIPT FOR NASA DATABASE 
% Author:
% Alberto Cuadra Lara, Universidad Carlos III de Madrid (UC3M)
% Last update: 12-Jun-2019 11:03

function data = read_CEA(filename)
% READ DATA FROM CEA AS TXT EXTENSION
% fid=fopen('test_soot_acetylene.txt','r'); 
fid=fopen(filename,'r'); 

ctLine=0;
i=0;
while ctLine<100000
    tline = fgetl(fid);
    if tline == -1, break, end
    if contains(tline,'CASE =')
        i = i+1;
        continue;
    end
    if contains(tline,'PHI,EQ.RATIO=')
        k = strfind(tline,'PHI,EQ.RATIO=');
        data.phi(i) = sscanf(tline(k+13:end),'%f');
    end
    if contains(tline,'P, BAR') 
        num = regexp(tline, '\d'); data.P(i) = sscanf(tline(num(1):num(end)),'%f'); tline = fgetl(fid);
        num = regexp(tline, '\d'); data.T(i) = sscanf(tline(num(1):num(end)),'%f'); tline = fgetl(fid);
        num = regexp(tline, '\d');
        if contains(tline,{'-'})
            idx = regexp(tline,'-');
            data.rho(i) = sscanf([tline(num(1):idx-1),'e',tline(idx:num(end))],'%f'); tline = fgetl(fid);
        else
            num = regexp(tline, '\d');
            num = sscanf(tline(num(1):num(end)),'%f');
            data.rho(i) = num(1); tline = fgetl(fid);
        end
        num = regexp(tline, '\d'); data.H(i) = sscanf(tline(num(1):num(end)),'%f'); tline = fgetl(fid);
        num = regexp(tline, '\d'); data.U(i) = sscanf(tline(num(1):num(end)),'%f'); tline = fgetl(fid);
        num = regexp(tline, '\d'); data.G(i) = sscanf(tline(num(1):num(end)),'%f'); tline = fgetl(fid);
        num = regexp(tline, '\d'); data.S(i) = sscanf(tline(num(1):num(end)),'%f'); tline = fgetl(fid);
        continue
    end
    if contains(tline,'Cp, KJ/(KG)(K)') 
        num = regexp(tline, '\d'); data.cp(i) = sscanf(tline(num(2):num(end)),'%f'); tline = fgetl(fid);
        continue
    end
    if contains(tline,'RHO/RHO1') 
        num = regexp(tline, '\d'); data.rho2rho1(i) = sscanf(tline(num(2):num(end)),'%f'); tline = fgetl(fid);
        tline = fgetl(fid);
        num = regexp(tline, '\d'); data.u1(i) = sscanf(tline(num(1):num(end)),'%f'); tline = fgetl(fid);
        continue
    end
    if contains(tline,'MOLE FRACTIONS')
        tline = fgetl(fid); tline = fgetl(fid);
        j=1;
        while ~contains(tline,'THERMODYNAMIC PROPERTIES FITTED TO 20000.K')
            if isempty(tline), break, end
            [sp1,sp2] = regexp(tline, '(?![*,-])\S\w*\s');  
            [mole,~] = regexp(tline, '\s\d');
            idx = regexp(tline,'-');
            if contains(tline,'C(gr)')
                data.X(i).mole{j,1} = 'Cbgrb';
            else
                data.X(i).mole{j,1} = tline(sp1:sp2-1);
            end
            data.X(i).mole{j,2} = sscanf([tline(mole:idx-1),'e',tline(idx:end)],'%f'); tline = fgetl(fid);
            j=j+1;
        end
    end
    ctLine=ctLine+1;
end
fclose(fid);