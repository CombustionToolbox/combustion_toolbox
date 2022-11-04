function gui_FLAG_IACValueChanged(app)
    if ~app.FLAG_IAC.Value
        app.text_P1.Text = 'FAC (only one option)';
        app.text_P1.Visible = 'on';
        app.text_RP1_2.Visible = 'on'; app.text_RP2_2.Visible = 'on';
        app.PP1.Visible = 'on'; app.PP2.Visible = 'on';
        app.PP1.Value = ''; app.PP2.Value = '';
    else
        app.text_P1.Text = 'Products';
        app.text_P1.Visible = 'off';
        app.text_RP1_2.Visible = 'off'; app.text_RP2_2.Visible = 'off';
        app.PP1.Visible = 'off'; app.PP2.Visible = 'off';
        app.PP1.Value = '2500'; app.PP2.Value = '1';
    end
end