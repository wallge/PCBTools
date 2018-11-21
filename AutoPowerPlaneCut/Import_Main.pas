
procedure TPlaneCut_Form.BtnFileClick(Sender: TObject);
begin
    If OpenDlg.Execute Then Begin
        TxtFile.Text := OpenDlg.FileName;
    End;
end;

procedure TPlaneCut_Form.BtnImportClick(Sender: TObject);
begin
    ImportPlaneCuts(TxtFile.Text);
    Close;
end;

