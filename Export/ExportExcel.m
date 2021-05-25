function ExportExcel(app)
    % delete(filename);
    if app.Misc.save_Excel
        ExcelOutputMatrix_R = FormattedOutput_test([],app.PD.phi.value,app.PS.strR,app.S.namespecies);
        ExcelOutputMatrix_P = FormattedOutput_test([],app.PD.phi.value,app.PS.strP,app.S.namespecies);
        [STATUS1,MESSAGE1]  = xlswrite(strcat(app.C.filename,'_R.xls'),ExcelOutputMatrix_R);
        [STATUS2,MESSAGE2]  = xlswrite(strcat(app.C.filename,'_P.xls'),ExcelOutputMatrix_P);
    end
    % [NUM,TXT,RAW] = xlsread(filename);
end