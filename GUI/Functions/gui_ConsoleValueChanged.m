function gui_ConsoleValueChanged(self, event)
    if strcmp('clear', self.Console.Value{1,1})
        ClearButtonPushed(self, event);
        output = ' ';
    else
        try
            output = evalc(self.Console.Value{1,1});
        catch error
            output = error.message;
        end
    end
    if ~isempty(output) || strcmp('clc', self.Console.Value{1,1})
        self.Console_text.Value = output;
    end
    self.Console.Value = '';
end