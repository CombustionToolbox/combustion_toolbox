function closing(strP,phi,display_species,timer_0,NameSpecies,mintol_display,ProblemType)
% fclose(output);
toc(timer_0);
% waitbar(1,f,'Finishing'); 
% answer = questdlg('Would you like to open the results?', ...
%     'Open results', ...
%     'Yes','No','Yes');
% Handle response
% switch answer
%     case 'Yes'
%         open('output.txt');
        if numel(phi)>1 && all(phi(2:end) ~= phi(1)) && ~strcmp(ProblemType,'DET_OVERDRIVEN')
            if isempty(display_species) 
                displaysweepresults(strP,phi,NameSpecies,mintol_display);
            else
                displaysweepresults(strP,phi,NameSpecies,mintol_display,display_species);
            end
        end
% end
% delete(f);