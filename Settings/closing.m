function closing(app,strP,phi,display_species,timer_0,LS,mintol_display,ProblemType)
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
                displaysweepresults(strP,phi,LS,mintol_display,'Equivalence Ratio $\phi$');
            else
                displaysweepresults(strP,phi,LS,mintol_display,'Equivalence Ratio $\phi$',display_species);
            end
            if ~any(strcmp(ProblemType,{'TP','TV'}))
                app.Misc.config.tit = ProblemType;
                app.Misc.config.labelx = 'Equivalence Ratio $\phi [-]$';
                app.Misc.config.labely = 'Temperature $T [K]$';
                plot_figure(app.PD.phi.value,app.PS.strP,'phi','T',app.Misc.config,app.PD.CompleteOrIncomplete);
            end
        elseif numel(phi)>1 && all(phi(2:end) == phi(1))
            app.Misc.config.tit = ProblemType;
            app.Misc.config.labelx = 'Critical Equivalence Ratio $\phi_c$';
            app.Misc.config.labely = 'Temperature $T [K]$';
            plot_figure(app.PS.strP,app.PS.strP,'phi_c','T',app.Misc.config,app.PD.CompleteOrIncomplete);
        end
% end
% delete(f);