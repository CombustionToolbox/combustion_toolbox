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
FLAG_ROW = false;

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
    
    if contains(tline,'INITIAL GAS')
        tline = fgetl(fid);
        num = regexp(tline, '\d'); data.M1(i, :) = sscanf(tline(num(2):num(end)),'%f'); tline = fgetl(fid);
        num = regexp(tline, '\d'); data.u1(i, :) = sscanf(tline(num(2):num(end)),'%f'); tline = fgetl(fid);
        num = regexp(tline, '\d'); data.P1(i, :) = sscanf(tline(num(1):num(end)),'%f'); tline = fgetl(fid);
        num = regexp(tline, '\d'); data.T1(i, :) = sscanf(tline(num(1):num(end)),'%f'); tline = fgetl(fid);
        
        if length(data.T1(i, :)) > 1
            FLAG_ROW = true;
        end

        num = regexp(tline, '\d');
        if contains(tline,{'-'})
            idx = regexp(tline,'-');
            for j = length(data.T1(i, :)):-1:1
                data.rho1(i, j) = sscanf([tline(num(j-1):idx(j-1)-1),'e',tline(idx(j):num(j))],'%f');
            end
        else
            num = regexp(tline, '\d');
            num = sscanf(tline(num(1):num(end)),'%f');
            num = num(num>0);
            data.rho1(i, :) = num;
        end
        tline = fgetl(fid);
        num = regexp(tline, '\d'); data.H1(i, :) = sscanf(tline(num(1)-1:num(end)),'%f'); tline = fgetl(fid);
        num = regexp(tline, '\d'); data.U1(i, :) = sscanf(tline(num(1)-1:num(end)),'%f'); tline = fgetl(fid);
        num = regexp(tline, '\d'); data.G1(i, :) = sscanf(tline(num(1)-1:num(end)),'%f'); tline = fgetl(fid);
        num = regexp(tline, '\d'); data.S1(i, :) = sscanf(tline(num(1):num(end)),'%f'); tline = fgetl(fid);
        tline = fgetl(fid);
        num = regexp(tline, '\d'); data.W1(i, :) = sscanf(tline(num(2):num(end)),'%f'); tline = fgetl(fid);
        num = regexp(tline, '\d'); data.cp1(i, :) = sscanf(tline(num(1):num(end)),'%f'); tline = fgetl(fid);
        num = regexp(tline, '\d'); data.gamma_s1(i, :) = sscanf(tline(num(1):num(end)),'%f'); tline = fgetl(fid);
        num = regexp(tline, '\d'); data.sound1(i, :) = sscanf(tline(num(1):num(end)),'%f'); tline = fgetl(fid);
        tline = fgetl(fid);
    end
    
    if contains(tline,'SHOCKED GAS (2)--INCIDENT')
        tline = fgetl(fid);
        tline = fgetl(fid);
        num = regexp(tline, '\d'); data.P2(i, :) = sscanf(tline(num(1):num(end)),'%f'); tline = fgetl(fid);
        num = regexp(tline, '\d'); data.T2(i, :) = sscanf(tline(num(1):num(end)),'%f'); tline = fgetl(fid);
        tline = fgetl(fid);

        data.rho2(i, :) = zeros(1, length(data.T2(i, :)));

        num = regexp(tline, '\d'); data.H2(i, :) = sscanf(tline(num(1)-1:num(end)),'%f'); tline = fgetl(fid);
        num = regexp(tline, '\d'); data.U2(i, :) = sscanf(tline(num(1)-1:num(end)),'%f'); tline = fgetl(fid);
        num = regexp(tline, '\d'); data.G2(i, :) = sscanf(tline(num(1)-1:num(end)),'%f'); tline = fgetl(fid);
        num = regexp(tline, '\d'); data.S2(i, :) = sscanf(tline(num(1):num(end)),'%f'); tline = fgetl(fid);
        tline = fgetl(fid);
        num = regexp(tline, '\d'); data.W2(i, :) = sscanf(tline(num(2):num(end)),'%f'); tline = fgetl(fid);
        num = regexp(tline, '\d'); data.dVdp_T(i, :) = sscanf(tline(num(2)-3:num(end)),'%f'); tline = fgetl(fid);
        num = regexp(tline, '\d'); data.dVdT_p(i, :) = sscanf(tline(num(2)-3:num(end)),'%f'); tline = fgetl(fid);
        num = regexp(tline, '\d'); data.cp2(i, :) = sscanf(tline(num(1):num(end)),'%f'); tline = fgetl(fid);
        num = regexp(tline, '\d'); data.gamma_s2(i, :) = sscanf(tline(num(1):num(end)),'%f'); tline = fgetl(fid);
        num = regexp(tline, '\d'); data.sound2(i, :) = sscanf(tline(num(1):num(end)),'%f'); tline = fgetl(fid);

        tline = fgetl(fid); tline = fgetl(fid); tline = fgetl(fid); tline = fgetl(fid); 
        num = regexp(tline, '\d'); data.rho2(i, :) = data.rho1(i, :) .* sscanf(tline(num(3)-1:num(end)),'%f')'; tline = fgetl(fid);
        num = regexp(tline, '\d'); data.u2(i, :) = sscanf(tline(num(2):num(end)),'%f'); tline = fgetl(fid);
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
        num = regexp(tline, '\d'); data.H(i) = sscanf(tline(num(1)-1:num(end)),'%f'); tline = fgetl(fid);
        num = regexp(tline, '\d'); data.U(i) = sscanf(tline(num(1)-1:num(end)),'%f'); tline = fgetl(fid);
        num = regexp(tline, '\d'); data.G(i) = sscanf(tline(num(1)-1:num(end)),'%f'); tline = fgetl(fid);
        num = regexp(tline, '\d'); data.S(i) = sscanf(tline(num(1):num(end)),'%f'); tline = fgetl(fid);
        continue
    end
    if contains(tline,'(dLV/dLP)t') 
        num = regexp(tline, '\d'); data.dVdp_T(i) = sscanf(tline(num(2)-2:num(end)),'%f'); tline = fgetl(fid);
        num = regexp(tline, '\d'); data.dVdT_p(i) = sscanf(tline(num(2)-2:num(end)),'%f'); tline = fgetl(fid);
    end

    if contains(tline,'Cp, KJ/(KG)(K)') 
        num = regexp(tline, '\d'); data.cp(i) = sscanf(tline(num(1):num(end)),'%f'); tline = fgetl(fid);
        if contains(tline,'GAMMAs') 
            num = regexp(tline, '\d'); data.gamma_s(i) = sscanf(tline(num(1):num(end)),'%f'); tline = fgetl(fid);
        end
        continue
    end
    
    if contains(tline,'RHO/RHO1') 
        num = regexp(tline, '\d'); data.rho2rho1(i) = sscanf(tline(num(3):num(end)),'%f'); tline = fgetl(fid);
        tline = fgetl(fid);
        num = regexp(tline, '\d'); data.u1(i) = sscanf(tline(num(1):num(end)),'%f'); tline = fgetl(fid);
        continue
    end

    if contains(tline,'MOLE FRACTIONS') && FLAG_ROW
        tline = fgetl(fid); tline = fgetl(fid);
        j=1;
        while ~contains(tline,'THERMODYNAMIC PROPERTIES FITTED TO 20000.K')
            if isempty(tline), break, end
            [sp1, sp2] = regexp(FullName2name(tline), '(?![*,-])\S\w*\s');
            [mole, ~] = regexp(tline, '\s\d');
            idx = regexp(tline,'-');
            
            for k = length(data.T2(i, :)):-1:1
                if contains(tline,'C(gr)')
                    data.X(k, 1).mole{j, 1} = 'Cbgrb';
                else
                    data.X(k, 1).mole{j, 1} = FullName2name(tline(sp1:sp2-1));
                end
                try
                    data.X(k, 1).mole{j, 2} = sscanf([tline(mole(k-1):idx(k-1)-1),'e',tline(idx(k-1):mole(k))],'%f');
                catch
                    % Last value
                    data.X(k, 1).mole{j, 2} = sscanf([tline(mole(k):idx(k)-1),'e',tline(idx(k):mole(k+1))],'%f');
                end
            end
            tline = fgetl(fid);
            j=j+1;
        end
    end

    if contains(tline,'MOLE FRACTIONS')
        tline = fgetl(fid); tline = fgetl(fid);
        j=1;
        while ~contains(tline,'THERMODYNAMIC PROPERTIES FITTED TO 20000.K')
            if isempty(tline), break, end
            [sp1, sp2] = regexp(FullName2name(tline), '(?![*,-])\S\w*\s');
            [mole, ~] = regexp(tline, '\s\d');
            idx = regexp(tline,'-');
            if contains(tline,'C(gr)')
                data.X(i, :).mole{j,1} = 'Cbgrb';
            else
                data.X(i, :).mole{j,1} = FullName2name(tline(sp1:sp2-1));
            end
            data.X(i, :).mole{j,2} = sscanf([tline(mole:idx-1),'e',tline(idx:end)],'%f'); tline = fgetl(fid);
            j=j+1;
        end
    end
    ctLine=ctLine+1;
end
fclose(fid);
end