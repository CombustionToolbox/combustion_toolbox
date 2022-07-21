function gui_FLAG_IACValueChanged(obj)
    if ~obj.FLAG_IAC.Value
        obj.text_P1.Text = 'FAC (only one option)';
        obj.text_P1.Visible = 'on';
        obj.text_RP1_2.Visible = 'on'; obj.text_RP2_2.Visible = 'on';
        obj.PP1.Visible = 'on'; obj.PP2.Visible = 'on';
        obj.PP1.Value = ''; obj.PP2.Value = '';
    else
        obj.text_P1.Text = 'Products';
        obj.text_P1.Visible = 'off';
        obj.text_RP1_2.Visible = 'off'; obj.text_RP2_2.Visible = 'off';
        obj.PP1.Visible = 'off'; obj.PP2.Visible = 'off';
        obj.PP1.Value = '2500'; obj.PP2.Value = '1';
    end
end